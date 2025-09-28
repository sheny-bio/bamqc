// 使用新IO适配子的main.rs示例片段

use clap::Parser;
use bamqc_io::{BamReader, BamError}; // 替换rust_htslib::bam
use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;
use tracing::{error, info, warn, debug};

/// 自定义错误类型（更新后）
#[derive(Error, Debug)]
enum BamQcError {
    #[error("BAM文件读取错误: {0}")]
    BamReadError(#[from] BamError), // 使用新的BamError
    #[error("过滤后没有可用于计算的配对读（TLEN>0的左端记录为空）")]
    NoValidReads,
    // ... 其他错误类型保持不变
}

/// 计算插入片段大小（使用新IO适配子）
fn compute_insert_size_with_new_io(
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

    // 使用新的统一API
    let mut reader = BamReader::from_path(bam_path)?;
    let mut stats = InsertSizeStats::new();

    info!("开始处理BAM文件: {}", bam_path);
    
    let mut processed_records = 0;
    let mut filtered_records = 0;

    // 使用新的迭代器接口，自动处理进度日志
    for record_result in reader.records() {
        let record = record_result?;
        processed_records += 1;

        // 所有的API调用保持不变，但现在使用封装的BamRecord
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
        
        if tlen <= 0 {
            continue;
        }

        let insert_size = tlen.abs() as i32;
        if insert_size == 0 {
            continue;
        }

        let orientation = determine_pair_orientation_new(&record);
        stats.add_insert_size(orientation, insert_size);
        filtered_records += 1;
    }

    // ... 后续处理逻辑保持不变
    Ok(0) // 示例返回值
}

/// 更新的配对方向判断函数
fn determine_pair_orientation_new(record: &bamqc_io::BamRecord) -> PairOrientation {
    let left_reverse = record.is_reverse();
    let right_reverse = record.is_mate_reverse();

    if left_reverse == right_reverse {
        PairOrientation::Tandem
    } else if !left_reverse && right_reverse {
        PairOrientation::Fr
    } else {
        PairOrientation::Rf
    }
}

// 其他结构体和枚举定义保持不变...
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, clap::ValueEnum)]
enum PairOrientation {
    #[clap(name = "fr")] Fr,
    #[clap(name = "rf")] Rf,
    #[clap(name = "tandem")] Tandem,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
enum Strategy { Specific, Dominant }

struct InsertSizeStats {
    histograms: HashMap<PairOrientation, HashMap<i32, u32>>,
    total_left_records: u32,
}

impl InsertSizeStats {
    fn new() -> Self {
        let mut histograms = HashMap::new();
        histograms.insert(PairOrientation::Fr, HashMap::new());
        histograms.insert(PairOrientation::Rf, HashMap::new());
        histograms.insert(PairOrientation::Tandem, HashMap::new());
        
        Self { histograms, total_left_records: 0 }
    }

    fn add_insert_size(&mut self, orientation: PairOrientation, size: i32) {
        *self.histograms.get_mut(&orientation).unwrap().entry(size).or_insert(0) += 1;
        self.total_left_records += 1;
    }
}
