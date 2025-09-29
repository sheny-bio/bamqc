use noodles::bam::{self, io::Reader};
use noodles::bgzf;
use noodles::sam::{self, alignment::record::Flags};
use std::fs::File;
use std::path::Path;
use thiserror::Error;
use tracing::info;

/// BAM/CRAM文件读取错误
#[derive(Error, Debug)]
pub enum BamError {
    #[error("IO错误: {0}")]
    IoError(#[from] std::io::Error),
    
    #[error("BAM格式错误: {0}")]
    BamError(String),
    
    #[error("SAM格式错误: {0}")]
    SamError(#[from] noodles::sam::header::ParseError),
    
    #[error("文件不存在: {path}")]
    FileNotFound { path: String },
}

/// BAM/CRAM文件读取器
pub struct BamReader {
    reader: Reader<bgzf::Reader<File>>,
    header: sam::Header,
    path: String,
}

impl std::fmt::Debug for BamReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BamReader")
            .field("path", &self.path)
            .field("header", &self.header)
            .finish()
    }
}

impl BamReader {
    /// 从文件路径创建BAM读取器
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, BamError> {
        let path_str = path.as_ref().to_string_lossy().to_string();
        
        if !path.as_ref().exists() {
            return Err(BamError::FileNotFound { path: path_str });
        }

        let file = File::open(&path)?;
        let mut reader = Reader::new(file);
        
        // 读取头部信息
        let header = reader.read_header().map_err(|e| BamError::BamError(e.to_string()))?;
        
        info!("已打开BAM文件: {}", path_str);
        
        Ok(Self {
            reader,
            header,
            path: path_str,
        })
    }

    /// 获取文件路径
    pub fn path(&self) -> &str {
        &self.path
    }

    /// 获取SAM头部信息
    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// 迭代所有记录
    pub fn records(&mut self) -> BamRecordIterator<'_> {
        BamRecordIterator {
            reader: &mut self.reader,
            count: 0,
        }
    }
}

/// BAM记录迭代器
pub struct BamRecordIterator<'a> {
    reader: &'a mut Reader<bgzf::Reader<File>>,
    count: u64,
}

impl<'a> Iterator for BamRecordIterator<'a> {
    type Item = Result<BamRecord, BamError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = bam::Record::default();
        match self.reader.read_record(&mut record) {
            Ok(0) => None, // EOF
            Ok(_) => {
                self.count += 1;
                Some(Ok(BamRecord { inner: record }))
            }
            Err(e) => Some(Err(BamError::BamError(e.to_string()))),
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
        self.inner.flags().contains(Flags::SEGMENTED)
    }

    /// 检查是否为secondary alignment
    #[inline]
    pub fn is_secondary(&self) -> bool {
        self.inner.flags().contains(Flags::SECONDARY)
    }

    /// 检查是否为supplementary alignment
    #[inline]
    pub fn is_supplementary(&self) -> bool {
        self.inner.flags().contains(Flags::SUPPLEMENTARY)
    }

    /// 检查是否标记为duplicate
    #[inline]
    pub fn is_duplicate(&self) -> bool {
        self.inner.flags().contains(Flags::DUPLICATE)
    }

    /// 检查是否未映射
    #[inline]
    pub fn is_unmapped(&self) -> bool {
        self.inner.flags().contains(Flags::UNMAPPED)
    }

    /// 检查mate是否未映射
    #[inline]
    pub fn is_mate_unmapped(&self) -> bool {
        self.inner.flags().contains(Flags::MATE_UNMAPPED)
    }

    /// 检查是否为proper pair
    #[inline]
    pub fn is_proper_pair(&self) -> bool {
        self.inner.flags().contains(Flags::PROPERLY_SEGMENTED)
    }

    /// 检查是否为reverse链
    #[inline]
    pub fn is_reverse(&self) -> bool {
        self.inner.flags().contains(Flags::REVERSE_COMPLEMENTED)
    }

    /// 检查mate是否为reverse链
    #[inline]
    pub fn is_mate_reverse(&self) -> bool {
        self.inner.flags().contains(Flags::MATE_REVERSE_COMPLEMENTED)
    }

    /// 获取参考序列ID
    #[inline]
    pub fn tid(&self) -> i32 {
        match self.inner.reference_sequence_id() {
            Some(Ok(id)) => id as i32,
            _ => -1,
        }
    }

    /// 获取mate的参考序列ID
    #[inline]
    pub fn mtid(&self) -> i32 {
        match self.inner.mate_reference_sequence_id() {
            Some(Ok(id)) => id as i32,
            _ => -1,
        }
    }

    /// 获取模板长度(insert size)
    #[inline]
    pub fn insert_size(&self) -> i64 {
        self.inner.template_length() as i64
    }

    /// 获取底层noodles记录的引用
    pub fn inner(&self) -> &bam::Record {
        &self.inner
    }
}