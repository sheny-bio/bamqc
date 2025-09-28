use rust_htslib::bam;
use std::path::Path;
use thiserror::Error;
use tracing::{debug, info};

/// BAM/CRAM文件读取错误
#[derive(Error, Debug)]
pub enum BamError {
    #[error("文件读取错误: {0}")]
    IoError(#[from] rust_htslib::errors::Error),
    
    #[error("文件不存在: {path}")]
    FileNotFound { path: String },
}

/// BAM/CRAM文件读取器
pub struct BamReader {
    reader: bam::Reader,
    path: String,
}

impl BamReader {
    /// 从文件路径创建BAM读取器
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, BamError> {
        let path_str = path.as_ref().to_string_lossy().to_string();
        
        if !path.as_ref().exists() {
            return Err(BamError::FileNotFound { path: path_str });
        }

        let reader = bam::Reader::from_path(&path)?;
        info!("已打开BAM文件: {}", path_str);
        
        Ok(Self {
            reader,
            path: path_str,
        })
    }

    /// 获取文件路径
    pub fn path(&self) -> &str {
        &self.path
    }

    /// 迭代所有记录
    pub fn records(&mut self) -> BamRecordIterator {
        BamRecordIterator {
            inner: self.reader.records(),
            count: 0,
        }
    }
}

/// BAM记录迭代器
pub struct BamRecordIterator<'a> {
    inner: bam::Records<'a, bam::Reader>,
    count: u64,
}

impl<'a> Iterator for BamRecordIterator<'a> {
    type Item = Result<BamRecord, BamError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.next() {
            Some(Ok(record)) => {
                self.count += 1;
                if self.count % 1_000_000 == 0 {
                    debug!("已处理 {} 条记录", self.count);
                }
                Some(Ok(BamRecord { inner: record }))
            }
            Some(Err(e)) => Some(Err(BamError::IoError(e))),
            None => None,
        }
    }
}

/// BAM记录封装
#[derive(Debug)]
pub struct BamRecord {
    inner: bam::Record,
}

impl BamRecord {
    /// 检查是否为配对读
    #[inline]
    pub fn is_paired(&self) -> bool {
        self.inner.is_paired()
    }

    /// 检查是否为secondary alignment
    #[inline]
    pub fn is_secondary(&self) -> bool {
        self.inner.is_secondary()
    }

    /// 检查是否为supplementary alignment
    #[inline]
    pub fn is_supplementary(&self) -> bool {
        self.inner.is_supplementary()
    }

    /// 检查是否标记为duplicate
    #[inline]
    pub fn is_duplicate(&self) -> bool {
        self.inner.is_duplicate()
    }

    /// 检查是否未映射
    #[inline]
    pub fn is_unmapped(&self) -> bool {
        self.inner.is_unmapped()
    }

    /// 检查mate是否未映射
    #[inline]
    pub fn is_mate_unmapped(&self) -> bool {
        self.inner.is_mate_unmapped()
    }

    /// 检查是否为proper pair
    #[inline]
    pub fn is_proper_pair(&self) -> bool {
        self.inner.is_proper_pair()
    }

    /// 检查是否为reverse链
    #[inline]
    pub fn is_reverse(&self) -> bool {
        self.inner.is_reverse()
    }

    /// 检查mate是否为reverse链
    #[inline]
    pub fn is_mate_reverse(&self) -> bool {
        self.inner.is_mate_reverse()
    }

    /// 获取参考序列ID
    #[inline]
    pub fn tid(&self) -> i32 {
        self.inner.tid()
    }

    /// 获取mate的参考序列ID
    #[inline]
    pub fn mtid(&self) -> i32 {
        self.inner.mtid()
    }

    /// 获取模板长度(insert size)
    #[inline]
    pub fn insert_size(&self) -> i64 {
        self.inner.insert_size()
    }

    /// 获取底层rust-htslib记录的引用
    pub fn inner(&self) -> &bam::Record {
        &self.inner
    }
}
