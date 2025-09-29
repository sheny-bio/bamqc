use clap::Parser;
use bamqc_io::bam::{BamReader, BamError};
use bamqc_core::{
    PairOrientation, Strategy, InsertSizeError, InsertSizeStats, 
    InsertSizeCalculator, determine_pair_orientation
};
use std::path::Path;
use tracing::{error, info, warn, debug};

/// 与Picard CollectInsertSizeMetrics一致的插入片段长度计算工具
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// 输入BAM/CRAM文件路径
    #[arg(short, long)]
    input: String,

    /// 包含标记为duplicate的读对
    #[arg(long)]
    include_duplicates: bool,

    /// 只统计proper pair
    #[arg(long)]
    require_proper_pair: bool,

    /// 类别最小占比阈值，丢弃占比低于该阈值的FR/RF/TANDEM类别
    #[arg(short = 'M', long, default_value = "0.05")]
    min_pct: f64,

    /// 配对方向类别
    #[arg(long, value_enum, default_value = "fr")]
    pair_orientation: PairOrientation,

    /// 输出策略
    #[arg(long, value_enum, default_value = "specific")]
    strategy: Strategy,

    /// 启用详细日志
    #[arg(short, long)]
    verbose: bool,
}


/// 自定义错误类型
#[derive(thiserror::Error, Debug)]
enum BamQcError {
    #[error("BAM文件读取错误: {0}")]
    BamReadError(#[from] BamError),
    #[error("插入片段大小计算错误: {0}")]
    InsertSizeError(#[from] InsertSizeError),
}




/// 计算插入片段大小
fn compute_insert_size(
    bam_path: &str,
    include_duplicates: bool,
    require_proper_pair: bool,
    min_pct: f64,
    orientation_pref: PairOrientation,
    strategy: Strategy,
) -> Result<i32, BamQcError> {
    let mut reader = BamReader::from_path(bam_path)?;
    let mut stats = InsertSizeStats::new();

    info!("开始处理BAM文件: {}", bam_path);
    
    let mut processed_records = 0;
    let mut filtered_records = 0;

    for result in reader.records() {
        let record = result?;
        processed_records += 1;

        if processed_records % 1_000_000 == 0 {
            debug!("已处理 {} 条记录", processed_records);
        }

        // 基础过滤
        if !record.is_paired() {
            continue;
        }
        if record.is_secondary() || record.is_supplementary() {
            continue;
        }
        if !include_duplicates && record.is_duplicate() {
            continue;
        }
        if record.is_unmapped() || record.is_mate_unmapped() {
            continue;
        }
        if record.tid() != record.mtid() {
            continue;
        }
        if require_proper_pair && !record.is_proper_pair() {
            continue;
        }

        let tlen = record.insert_size();
        
        // 只计"左端记录"（TLEN > 0）
        if tlen <= 0 {
            continue;
        }

        let insert_size = tlen.abs() as i32;
        if insert_size == 0 {
            continue;
        }

        let orientation = determine_pair_orientation(record.is_reverse(), record.is_mate_reverse());
        stats.add_insert_size(orientation, insert_size);
        filtered_records += 1;
    }

    info!("处理完成：总记录数 {}，有效左端记录数 {}", processed_records, filtered_records);

    // 使用 InsertSizeCalculator 来计算最终结果
    let result = InsertSizeCalculator::calculate(&stats, min_pct, orientation_pref, strategy)?;
    
    // 记录保留的类别信息
    for (orientation, counts) in &stats.histograms {
        let count: u32 = counts.values().sum();
        if count == 0 {
            continue;
        }
        let pct = count as f64 / stats.total_left_records as f64;
        if pct >= min_pct {
            info!("保留类别 {}: {} 个读对 ({:.2}%)", orientation, count, pct * 100.0);
        } else {
            warn!("丢弃类别 {}: {} 个读对 ({:.2}%) < {:.1}%", orientation, count, pct * 100.0, min_pct * 100.0);
        }
    }

    match strategy {
        Strategy::Specific => {
            info!("使用指定方向 {} 的中位数: {}", orientation_pref, result);
        }
        Strategy::Dominant => {
            // 找到最大类别
            let best_orientation = stats.histograms
                .iter()
                .filter_map(|(orientation, counts)| {
                    let count: u32 = counts.values().sum();
                    let pct = count as f64 / stats.total_left_records as f64;
                    if pct >= min_pct {
                        Some((orientation, count))
                    } else {
                        None
                    }
                })
                .max_by_key(|(_, count)| *count)
                .map(|(orientation, _)| orientation);
            
            if let Some(best_orientation) = best_orientation {
                info!("使用最大类别 {} 的中位数: {}", best_orientation, result);
            }
        }
    }

    Ok(result)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // 初始化日志
    let log_level = if args.verbose { "debug" } else { "info" };
    tracing_subscriber::fmt()
        .with_env_filter(format!("bamqc={}", log_level))
        .init();

    // 验证输入文件存在
    if !Path::new(&args.input).exists() {
        error!("输入文件不存在: {}", args.input);
        std::process::exit(1);
    }

    match compute_insert_size(
        &args.input,
        args.include_duplicates,
        args.require_proper_pair,
        args.min_pct,
        args.pair_orientation,
        args.strategy,
    ) {
        Ok(median_size) => {
            // 按需求：只输出一个整数
            println!("{}", median_size);
            Ok(())
        }
        Err(e) => {
            error!("{}", e);
            std::process::exit(1);
        }
    }
}
