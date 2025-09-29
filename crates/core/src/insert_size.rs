//! BAM质量控制的插入片段大小计算模块。
//! 
//! 本模块提供从配对末端测序数据计算插入片段大小的功能，
//! 支持不同的配对方向和计算策略。

use std::collections::HashMap;
use thiserror::Error;

/// 插入片段大小计算的配对方向类型。
/// 
/// 表示配对末端读长在参考基因组中相对于彼此的不同方向。
/// 
/// # Examples
/// 
/// ```rust
/// use bamqc_core::PairOrientation;
/// 
/// let orientation = PairOrientation::Fr;
/// assert_eq!(orientation.to_string(), "FR");
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, clap::ValueEnum)]
pub enum PairOrientation {
    /// Forward-Reverse方向（典型的文库制备方式）。
    /// 
    /// 第一个读长比对到正向链，第二个读长比对到反向链。
    #[clap(name = "fr")]
    Fr,
    /// Reverse-Forward方向。
    /// 
    /// 第一个读长比对到反向链，第二个读长比对到正向链。
    #[clap(name = "rf")]
    Rf,
    /// 串联方向。
    /// 
    /// 两个读长都比对到同一条链（要么都是正向，要么都是反向）。
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
/// 
/// 当数据中存在多种配对方向时，此枚举决定使用哪一种
/// 来进行最终的插入片段大小计算。
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum Strategy {
    /// 使用指定的配对方向类别。
    /// 
    /// 如果指定的方向不满足最小百分比阈值，将返回错误。
    Specific,
    /// 使用主导的配对方向类别。
    /// 
    /// 在满足最小百分比阈值的方向中选择读长数量最多的那个。
    Dominant,
}

/// 插入片段大小计算过程中可能发生的错误。
/// 
/// 这些错误表示从配对末端测序数据计算插入片段大小时的
/// 各种失败模式。
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
}

/// 插入片段大小统计结果。
/// 
/// 存储不同配对方向的插入片段大小分布统计信息。
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


pub fn calculate_insert_size(