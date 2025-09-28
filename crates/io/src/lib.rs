//! BAM/CRAM文件IO适配子库
//!
//! 提供统一的BAM/CRAM文件读取API，封装rust-htslib的复杂性

pub mod bam;

// 重新导出主要类型
pub use bam::{BamError, BamReader, BamRecord, BamRecordIterator};

/// 库版本信息
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// 预加载函数，便于快速创建BamReader
pub fn open_bam<P: AsRef<std::path::Path>>(path: P) -> Result<BamReader, BamError> {
    BamReader::from_path(path)
}
