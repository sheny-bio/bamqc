#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
picard_insert_size_cli.py

与 Picard CollectInsertSizeMetrics 口径一致地计算插入片段长度（MEDIAN_INSERT_SIZE），仅打印一个整数。

口径要点（与 Picard 对齐）：
1) 仅使用 SAM TLEN 字段：insert_size = abs(TLEN)
2) 只统计每对里的“左端记录”：即 TLEN > 0 的那条；避免双倍计数
3) 过滤逻辑：
   - 仅成对 (is_paired)
   - 双端均已比对 (not is_unmapped, not mate_is_unmapped)
   - 同参考序列 (reference_id == next_reference_id)
   - 排除二级/补充比对 (not is_secondary, not is_supplementary)
   - 默认排除 duplicates（可用 --include-duplicates 改为包含）
   - 默认不强制 proper pair（Picard 通常不硬性要求；可用 --require-proper-pair 开启）
4) 按 FR/RF/TANDEM 三类分桶，并按 MINIMUM_PCT 丢弃占比低的类别（默认 0.05）
5) 输出策略：
   - 默认：输出 FR 类别（若被丢弃则报错）
   - 或者选择：--pair-orientation RF/TANDEM
   - 或者使用：--strategy dominant  输出保留类别中“读数最多的那一类”的中位数

仅打印一个整数，便于下游脚本对接。
"""

import argparse
import sys
from collections import Counter
import pysam


def _pair_orientation_leftmost(rec: pysam.AlignedSegment) -> str:
    """
    仅在 TLEN > 0（该记录为左端）时调用。
    Picard/HTSJDK 的 FR/RF/TANDEM 语义可用“左端/右端的链方向”来表达：
      - 左端为正(+) 且右端为负(-) -> FR（常见文库）
      - 左端为负(-) 且右端为正(+) -> RF
      - 同向(++, --) -> TANDEM
    这里：rec 是左端；mate 是右端。
    """
    left_rev = rec.is_reverse             # 左端（当前记录）
    right_rev = rec.mate_is_reverse       # 右端（mate）

    if left_rev == right_rev:
        return "TANDEM"
    return "FR" if (not left_rev and right_rev) else "RF"


def _hist_median_from_counts(counts: Counter) -> int:
    """
    按直方图“累计频数首次 >= 50%”所在 bin 的 key 作为中位数（整数）。
    与 Picard/HTSJDK 的 Histogram 分位实现一致（不取两数均值）。
    """
    if not counts:
        return 0
    total = sum(counts.values())
    # “上中位”门槛：1-based 计数
    threshold = (total + 1) // 2
    running = 0
    for size in sorted(counts.keys()):
        running += counts[size]
        if running >= threshold:
            return int(size)
    # 理论不可达
    return int(max(counts.keys()))


def compute_insert_size(
    bam_path: str,
    include_duplicates: bool = False,
    require_proper_pair: bool = False,
    min_pct: float = 0.05,
    orientation_pref: str = "FR",
    strategy: str = "specific",
) -> int:
    """
    返回一个整数 insert_size（MEDIAN_INSERT_SIZE）。
    strategy:
      - "specific": 仅输出 orientation_pref（默认 FR）
      - "dominant": 输出保留类别中“读数最多的那一类”的中位数
    """
    if strategy not in {"specific", "dominant"}:
        raise ValueError("strategy 必须为 'specific' 或 'dominant'")
    if orientation_pref not in {"FR", "RF", "TANDEM"}:
        raise ValueError("pair-orientation 必须为 FR/RF/TANDEM")

    hists = {"FR": Counter(), "RF": Counter(), "TANDEM": Counter()}
    total_left_records = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            # 基础过滤
            if not rec.is_paired:
                continue
            if rec.is_secondary or rec.is_supplementary:
                continue
            if (not include_duplicates) and rec.is_duplicate:
                continue
            if rec.is_unmapped or rec.mate_is_unmapped:
                continue
            if rec.reference_id != rec.next_reference_id:
                continue
            if require_proper_pair and (not rec.is_proper_pair):
                continue

            tlen = rec.template_length
            # 只计“左端记录”（TLEN > 0）
            if tlen <= 0:
                continue

            ins = abs(int(tlen))
            if ins == 0:
                continue

            ori = _pair_orientation_leftmost(rec)
            hists[ori][ins] += 1
            total_left_records += 1

    if total_left_records == 0:
        raise ValueError("过滤后没有可用于计算的配对读（TLEN>0 的左端记录为空）。")

    # MINIMUM_PCT：按类别占比丢弃
    kept = {}
    for ori, cnts in hists.items():
        n = sum(cnts.values())
        if n == 0:
            continue
        pct = n / total_left_records
        if pct >= min_pct:
            kept[ori] = (n, cnts)

    if not kept:
        raise ValueError(
            f"所有方向类别占比均 < MINIMUM_PCT={min_pct}，无法给出 insert_size。"
        )

    if strategy == "specific":
        if orientation_pref not in kept:
            raise ValueError(
                f"所选方向 {orientation_pref} 被 MINIMUM_PCT={min_pct} 丢弃，"
                f"可降低阈值或改用 --strategy dominant。"
            )
        counts = kept[orientation_pref][1]
        return _hist_median_from_counts(counts)

    # strategy == "dominant"
    best_ori = max(kept.items(), key=lambda kv: kv[1][0])[0]
    counts = kept[best_ori][1]
    return _hist_median_from_counts(counts)


def main():
    ap = argparse.ArgumentParser(
        description="Picard CollectInsertSizeMetrics 口径的插入片段长度（MEDIAN_INSERT_SIZE）计算器；仅打印一个整数。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("-i", "--input", required=True, help="输入 BAM/CRAM（建议 BAM）")
    ap.add_argument(
        "--include-duplicates",
        action="store_true",
        help="包含标记为 duplicate 的读对（Picard 默认不包含）。",
    )
    ap.add_argument(
        "--require-proper-pair",
        action="store_true",
        help="只统计 proper pair（Picard 通常不强制；如需更严格可开启）。",
    )
    ap.add_argument(
        "-M", "--min-pct",
        type=float, default=0.05,
        help="类别最小占比阈值（MINIMUM_PCT），丢弃占比低于该阈值的 FR/RF/TANDEM 类别。",
    )
    ap.add_argument(
        "--pair-orientation",
        choices=["FR", "RF", "TANDEM"],
        default="FR",
        help="在 strategy=specific 下要输出的配对方向类别。",
    )
    ap.add_argument(
        "--strategy",
        choices=["specific", "dominant"],
        default="specific",
        help="specific：输出 --pair-orientation 指定类别；dominant：输出保留类别中读数最多的那一类。",
    )

    args = ap.parse_args()

    if not (0.0 <= args.min_pct <= 0.5):
        print("ERROR: --min-pct 必须在 [0, 0.5] 之间。", file=sys.stderr)
        sys.exit(2)

    try:
        median_size = compute_insert_size(
            bam_path=args.input,
            include_duplicates=args.include_duplicates,
            require_proper_pair=args.require_proper_pair,
            min_pct=args.min_pct,
            orientation_pref=args.pair_orientation,
            strategy=args.strategy,
        )
        # 按需求：只输出一个整数
        print(median_size)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
