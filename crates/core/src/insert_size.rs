//! BAM质量控制的插入片段大小计算模块。
//! 
//! 本模块提供从配对末端测序数据计算插入片段大小的功能，
//! 支持不同的配对方向和计算策略。

use std::collections::HashMap;
use thiserror::Error;
use bamqc_io::bam::{BamReader, BamError};
use tracing::{info, warn, debug};

/// 插入片段大小计算的配对方向类型。
/// 
/// 表示配对末端读长在参考基因组中相对于彼此的不同方向。
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, clap::ValueEnum)]
pub enum PairOrientation {
    /// Forward-Reverse方向（典型的文库制备方式）。
    #[clap(name = "fr")]
    Fr,
    /// Reverse-Forward方向。
    #[clap(name = "rf")]
    Rf,
    /// 串联方向。
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

/// 选择用于插入片段大小计算的配对方向策略。
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum Strategy {
    /// 使用指定的配对方向类别。
    Specific,
    /// 使用主导的配对方向类别。
    Dominant,
}

/// 插入片段大小计算过程中可能发生的错误。
#[derive(Error, Debug)]
pub enum InsertSizeError {
    /// 没有可用于计算的有效配对读长。
    /// 
    /// 当过滤后没有TLEN > 0的左端记录时发生。
    #[error("过滤后没有可用于计算的配对读（TLEN>0的左端记录为空）")]
    NoValidReads,
    
    /// 所有方向类别都被最小百分比阈值过滤掉了。
    /// 
    /// 当没有任何配对方向类别满足最小百分比要求时发生。
    #[error("所有方向类别占比均 < MINIMUM_PCT={min_pct:.3}，无法给出insert_size")]
    AllCategoriesFiltered { 
        /// 应用的最小百分比阈值
        min_pct: f64 
    },
    
    /// 指定的方向被最小百分比阈值过滤掉了。
    /// 
    /// 当使用`Specific`策略且请求的方向不满足最小百分比要求时发生。
    #[error("所选方向 {orientation} 被MINIMUM_PCT={min_pct:.3}丢弃，可降低阈值或改用dominant策略")]
    OrientationFiltered {
        /// 被过滤掉的方向
        orientation: PairOrientation,
        /// 应用的最小百分比阈值
        min_pct: f64,
    },
    
    /// 无效的最小百分比值。
    /// 
    /// 最小百分比必须在0.0和0.5之间（含边界值）。
    #[error("min_pct必须在[0, 0.5]之间")]
    InvalidMinPct,
    
    /// BAM文件IO错误。
    /// 
    /// 当读取BAM文件时发生IO错误时发生。
    #[error("BAM文件IO错误: {0}")]
    BamError(#[from] BamError),
}

/// 插入片段大小统计结果。
#[derive(Debug)]
pub struct InsertSizeStats {
    /// 各方向的插入大小计数直方图。
    pub histograms: HashMap<PairOrientation, HashMap<i32, u32>>,
    
    /// 总的左端记录数。
    pub total_left_records: u32,
}

impl InsertSizeStats {

    pub fn new() -> Self {
        let mut histograms = HashMap::new();
        histograms.insert(PairOrientation::Fr, HashMap::new());
        histograms.insert(PairOrientation::Rf, HashMap::new());
        histograms.insert(PairOrientation::Tandem, HashMap::new());
        
        Self {
            histograms,
            total_left_records: 0,
        }
    }

    /// 添加一个插入大小记录。
    /// 
    /// # Parameters
    /// 
    /// * `orientation` - 配对方向类型
    /// * `size` - 插入片段大小
    /// 
    /// # Notes
    /// 
    /// 同时更新对应方向的直方图计数和总记录数。
    pub fn add_insert_size(&mut self, orientation: PairOrientation, size: i32) {
        *self.histograms.get_mut(&orientation).unwrap().entry(size).or_insert(0) += 1;
        self.total_left_records += 1;
    }
}

impl Default for InsertSizeStats {
    fn default() -> Self {
        Self::new()
    }
}

/// 确定配对方向（仅在TLEN > 0时调用）。
/// 
/// 根据Picard/HTSJDK的FR/RF/TANDEM语义确定配对读长的方向类型。
/// 
/// # Parameters
/// 
/// * `left_reverse` - 左端读长是否为反向
/// * `right_reverse` - 右端读长是否为反向
/// 
/// # Returns
/// 
/// 返回对应的配对方向类型：
/// - 左端为正(+)且右端为负(-) -> FR（常见文库）
/// - 左端为负(-)且右端为正(+) -> RF  
/// - 同向(++, --) -> TANDEM
/// 
/// # Examples
/// 
/// ```
/// use bamqc_core::{determine_pair_orientation, PairOrientation};
/// 
/// // 典型的FR文库
/// assert_eq!(determine_pair_orientation(false, true), PairOrientation::Fr);
/// // RF文库
/// assert_eq!(determine_pair_orientation(true, false), PairOrientation::Rf);
/// // 串联方向
/// assert_eq!(determine_pair_orientation(false, false), PairOrientation::Tandem);
/// ```
pub fn determine_pair_orientation(left_reverse: bool, right_reverse: bool) -> PairOrientation {
    if left_reverse == right_reverse {
        PairOrientation::Tandem
    } else if !left_reverse && right_reverse {
        PairOrientation::Fr
    } else {
        PairOrientation::Rf
    }
}


/// 插入片段大小计算器。
/// 
/// 提供从统计数据计算插入片段大小的静态方法。
pub struct InsertSizeCalculator;

impl InsertSizeCalculator {
    /// 从计数HashMap计算中位数。
    /// 
    /// 按直方图“累计频数首次 >= 50%”所在bin的key作为中位数（整数）。
    /// 与Picard/HTSJDK的Histogram分位实现一致（不取两数均值）。
    /// 
    /// # Parameters
    /// 
    /// * `counts` - 插入大小到出现次数的映射
    /// 
    /// # Returns
    /// 
    /// 返回计算得到的中位数，如果输入为空则返回0。
    pub fn calculate_median_from_counts(counts: &HashMap<i32, u32>) -> i32 {
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
        
        // 理论上不可达，返回最大值作为备选
        *counts.keys().max().unwrap_or(&0)
    }

    /// 从统计数据计算最终的插入片段大小。
    /// 
    /// 根据指定的策略和最小百分比阈值，从统计数据中选择合适的配对方向
    /// 并计算其中位数作为最终的插入片段大小。
    /// 
    /// # Parameters
    /// 
    /// * `stats` - 插入大小统计数据
    /// * `min_pct` - 最小百分比阈值（必须在0.0-0.5之间）
    /// * `orientation_pref` - 首选的配对方向（在Specific策略下使用）
    /// * `strategy` - 选择策略（Specific或Dominant）
    /// 
    /// # Returns
    /// 
    /// 成功时返回计算得到的插入片段大小，失败时返回相应错误。
    /// 
    /// # Errors
    /// 
    /// * `InvalidMinPct` - 当min_pct不在有效范围内时
    /// * `NoValidReads` - 当没有有效记录时
    /// * `AllCategoriesFiltered` - 当所有方向都被过滤时
    /// * `OrientationFiltered` - 当指定方向被过滤时（Specific策略）
    pub fn calculate(
        stats: &InsertSizeStats,
        min_pct: f64,
        orientation_pref: PairOrientation,
        strategy: Strategy,
    ) -> Result<i32, InsertSizeError> {
        if !(0.0..=0.5).contains(&min_pct) {
            return Err(InsertSizeError::InvalidMinPct);
        }

        if stats.total_left_records == 0 {
            return Err(InsertSizeError::NoValidReads);
        }

        // 按最小百分比阈值过滤方向类别
        let mut kept_categories = HashMap::new();
        for (orientation, counts) in &stats.histograms {
            let count: u32 = counts.values().sum();
            if count == 0 {
                continue;
            }
            let pct = count as f64 / stats.total_left_records as f64;
            if pct >= min_pct {
                kept_categories.insert(*orientation, (count, counts));
            }
        }

        if kept_categories.is_empty() {
            return Err(InsertSizeError::AllCategoriesFiltered { min_pct });
        }

        match strategy {
            Strategy::Specific => {
                if let Some((_, counts)) = kept_categories.get(&orientation_pref) {
                    let median = Self::calculate_median_from_counts(counts);
                    Ok(median)
                } else {
                    Err(InsertSizeError::OrientationFiltered {
                        orientation: orientation_pref,
                        min_pct,
                    })
                }
            }
            Strategy::Dominant => {
                let (_, (_, counts)) = kept_categories
                    .iter()
                    .max_by_key(|(_, (count, _))| *count)
                    .unwrap();
                let median = Self::calculate_median_from_counts(counts);
                Ok(median)
            }
        }
    }
}

/// 计算插入片段大小。
/// 
/// 从BAM文件中读取配对末端测序数据，计算插入片段大小的中位数。
/// 支持多种过滤条件和计算策略。
/// 
/// # Parameters
/// 
/// * `bam_path` - BAM文件路径
/// * `include_duplicates` - 是否包含标记为duplicate的读对
/// * `require_proper_pair` - 是否只统计proper pair
/// * `min_pct` - 类别最小占比阈值，丢弃占比低于该阈值的FR/RF/TANDEM类别
/// * `orientation_pref` - 首选的配对方向
/// * `strategy` - 输出策略（Specific或Dominant）
/// 
/// # Returns
/// 
/// 成功时返回计算得到的插入片段大小中位数，失败时返回相应错误。
pub fn compute_insert_size(
    bam_path: &str,
    include_duplicates: bool,
    require_proper_pair: bool,
    min_pct: f64,
    orientation_pref: PairOrientation,
    strategy: Strategy,
) -> Result<i32, InsertSizeError> {
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
        if !record.is_segmented() {
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
        if require_proper_pair && !record.is_properly_segmented() {
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