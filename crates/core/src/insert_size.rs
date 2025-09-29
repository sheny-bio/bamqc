use std::collections::HashMap;
use thiserror::Error;

/// 配对方向枚举
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, clap::ValueEnum)]
pub enum PairOrientation {
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
pub enum Strategy {
    /// 输出指定的配对方向类别
    Specific,
    /// 输出保留类别中读数最多的那一类
    Dominant,
}

/// 插入片段大小计算相关错误
#[derive(Error, Debug)]
pub enum InsertSizeError {
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
pub struct InsertSizeStats {
    /// 各方向的插入大小计数
    pub histograms: HashMap<PairOrientation, HashMap<i32, u32>>,
    /// 总的左端记录数
    pub total_left_records: u32,
}

impl InsertSizeStats {
    /// 创建新的统计实例
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

    /// 添加一个插入大小记录
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

/// 确定配对方向（仅在TLEN > 0时调用）
/// 
/// 根据Picard/HTSJDK的FR/RF/TANDEM语义：
/// - 左端为正(+)且右端为负(-) -> FR（常见文库）
/// - 左端为负(-)且右端为正(+) -> RF  
/// - 同向(++, --) -> TANDEM
pub fn determine_pair_orientation(left_reverse: bool, right_reverse: bool) -> PairOrientation {
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
    
    // 理论不可达
    *counts.keys().max().unwrap_or(&0)
}

/// 插入片段大小计算器
pub struct InsertSizeCalculator;

impl InsertSizeCalculator {
    /// 从统计数据计算最终的插入片段大小
    pub fn calculate_from_stats(
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
            }
        }

        if kept_categories.is_empty() {
            return Err(InsertSizeError::AllCategoriesFiltered { min_pct });
        }

        match strategy {
            Strategy::Specific => {
                if let Some((_, counts)) = kept_categories.get(&orientation_pref) {
                    let median = calculate_median_from_counts(counts);
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
                let median = calculate_median_from_counts(counts);
                Ok(median)
            }
        }
    }
}