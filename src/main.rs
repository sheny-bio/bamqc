use clap::Parser;
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;
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

/// 配对方向枚举
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, clap::ValueEnum)]
enum PairOrientation {
    /// Forward-Reverse方向（常见文库）
    #[clap(name = "fr")]
    Fr,
    /// Reverse-Forward方向
    #[clap(name = "rf")]
    Rf,
    /// 同向配对
    #[clap(name = "tandem")]
    Tandem,
}

impl std::fmt::Display for PairOrientation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PairOrientation::Fr => write!(f, "FR"),
            PairOrientation::Rf => write!(f, "RF"),
            PairOrientation::Tandem => write!(f, "TANDEM"),
        }
    }
}

/// 输出策略枚举
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
enum Strategy {
    /// 输出指定的配对方向类别
    Specific,
    /// 输出保留类别中读数最多的那一类
    Dominant,
}

/// 自定义错误类型
#[derive(Error, Debug)]
enum BamQcError {
    #[error("BAM文件读取错误: {0}")]
    BamReadError(#[from] rust_htslib::errors::Error),
    #[error("过滤后没有可用于计算的配对读（TLEN>0的左端记录为空）")]
    NoValidReads,
    #[error("所有方向类别占比均 < MINIMUM_PCT={min_pct:.3}，无法给出insert_size")]
    AllCategoriesFiltered { min_pct: f64 },
    #[error("所选方向 {orientation} 被MINIMUM_PCT={min_pct:.3}丢弃，可降低阈值或改用dominant策略")]
    OrientationFiltered {
        orientation: PairOrientation,
        min_pct: f64,
    },
    #[error("min_pct必须在[0, 0.5]之间")]
    InvalidMinPct,
}

/// 插入大小统计结果
#[derive(Debug)]
struct InsertSizeStats {
    /// 各方向的插入大小计数
    histograms: HashMap<PairOrientation, HashMap<i32, u32>>,
    /// 总的左端记录数
    total_left_records: u32,
}

impl InsertSizeStats {
    /// 创建新的统计实例
    fn new() -> Self {
        let mut histograms = HashMap::new();
        histograms.insert(PairOrientation::Fr, HashMap::new());
        histograms.insert(PairOrientation::Rf, HashMap::new());
        histograms.insert(PairOrientation::Tandem, HashMap::new());
        
        Self {
            histograms,
            total_left_records: 0,
        }
    }

    /// 添加一个插入大小记录
    fn add_insert_size(&mut self, orientation: PairOrientation, size: i32) {
        *self.histograms.get_mut(&orientation).unwrap().entry(size).or_insert(0) += 1;
        self.total_left_records += 1;
    }
}

/// 确定配对方向（仅在TLEN > 0时调用）
/// 
/// 根据Picard/HTSJDK的FR/RF/TANDEM语义：
/// - 左端为正(+)且右端为负(-) -> FR（常见文库）
/// - 左端为负(-)且右端为正(+) -> RF  
/// - 同向(++, --) -> TANDEM
fn determine_pair_orientation(record: &bam::Record) -> PairOrientation {
    let left_reverse = record.is_reverse();     // 左端（当前记录）
    let right_reverse = record.is_mate_reverse(); // 右端（mate）

    if left_reverse == right_reverse {
        PairOrientation::Tandem
    } else if !left_reverse && right_reverse {
        PairOrientation::Fr
    } else {
        PairOrientation::Rf
    }
}

/// 从计数HashMap计算中位数
/// 
/// 按直方图"累计频数首次 >= 50%"所在bin的key作为中位数（整数）
/// 与Picard/HTSJDK的Histogram分位实现一致（不取两数均值）
fn calculate_median_from_counts(counts: &HashMap<i32, u32>) -> i32 {
    if counts.is_empty() {
        return 0;
    }

    let total: u32 = counts.values().sum();
    let threshold = (total + 1) / 2; // "上中位"门槛：1-based计数
    
    let mut sorted_sizes: Vec<i32> = counts.keys().copied().collect();
    sorted_sizes.sort();
    
    let mut running = 0;
    for size in sorted_sizes {
        running += counts[&size];
        if running >= threshold {
            return size;
        }
    }
    
    // 理论不可达
    *counts.keys().max().unwrap_or(&0)
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
    if !(0.0..=0.5).contains(&min_pct) {
        return Err(BamQcError::InvalidMinPct);
    }

    let mut reader = bam::Reader::from_path(bam_path)?;
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

        let orientation = determine_pair_orientation(&record);
        stats.add_insert_size(orientation, insert_size);
        filtered_records += 1;
    }

    info!("处理完成：总记录数 {}，有效左端记录数 {}", processed_records, filtered_records);

    if stats.total_left_records == 0 {
        return Err(BamQcError::NoValidReads);
    }

    // 按MINIMUM_PCT过滤类别
    let mut kept_categories = HashMap::new();
    for (orientation, counts) in &stats.histograms {
        let count: u32 = counts.values().sum();
        if count == 0 {
            continue;
        }
        let pct = count as f64 / stats.total_left_records as f64;
        if pct >= min_pct {
            kept_categories.insert(*orientation, (count, counts));
            info!("保留类别 {}: {} 个读对 ({:.2}%)", orientation, count, pct * 100.0);
        } else {
            warn!("丢弃类别 {}: {} 个读对 ({:.2}%) < {:.1}%", orientation, count, pct * 100.0, min_pct * 100.0);
        }
    }

    if kept_categories.is_empty() {
        return Err(BamQcError::AllCategoriesFiltered { min_pct });
    }

    match strategy {
        Strategy::Specific => {
            if let Some((_, counts)) = kept_categories.get(&orientation_pref) {
                let median = calculate_median_from_counts(counts);
                info!("使用指定方向 {} 的中位数: {}", orientation_pref, median);
                Ok(median)
            } else {
                Err(BamQcError::OrientationFiltered {
                    orientation: orientation_pref,
                    min_pct,
                })
            }
        }
        Strategy::Dominant => {
            let (best_orientation, (_, counts)) = kept_categories
                .iter()
                .max_by_key(|(_, (count, _))| *count)
                .unwrap();
            let median = calculate_median_from_counts(counts);
            info!("使用最大类别 {} 的中位数: {}", best_orientation, median);
            Ok(median)
        }
    }
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
