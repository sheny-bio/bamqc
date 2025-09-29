use noodles::bam::{self, io::Reader};
use noodles::bgzf;
use noodles::sam::{self};
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

    /// 是否为配对读; 对应flag: 0x1
    pub fn is_segmented(&self) -> bool {
        self.inner.flags().is_segmented()
    }

    
    /// 是否为主要读对; 对应flag: 0x2
    pub fn is_properly_segmented(&self) -> bool {
        self.inner.flags().is_properly_segmented()
    }

    /// reads未能比对到参考序列; 对应flag: 0x4
    pub fn is_unmapped(&self) -> bool {
        self.inner.flags().is_unmapped()
    }

    /// 配对reads的mate未能比对到参考序列; 对应flag: 0x8
    pub fn is_mate_unmapped(&self) -> bool {
        self.inner.flags().is_mate_unmapped()
    }

    /// reads是否为反向链; 对应flag: 0x10
    pub fn is_reverse_complemented(&self) -> bool {
        self.inner.flags().is_reverse_complemented()
    }

    pub fn is_reverse(&self) -> bool {
        self.inner.flags().is_reverse_complemented()
    }

    /// 配对reads的mate是否为反向链; 对应flag: 0x20
    pub fn is_mate_reverse_complemented(&self) -> bool {
        self.inner.flags().is_mate_reverse_complemented()
    }

    pub fn is_mate_reverse(&self) -> bool {
        self.inner.flags().is_mate_reverse_complemented()
    }
    
    /// reads是否在5'端; 对应flag: 0x40
    pub fn is_first_segment(&self) -> bool {
        self.inner.flags().is_first_segment()
    }

    /// reads是否在3'端; 对应flag: 0x80
    pub fn is_last_segment(&self) -> bool {
        self.inner.flags().is_last_segment()
    }

    /// 是否为次要比对; 对应flag: 0x100
    pub fn is_secondary(&self) -> bool {
        self.inner.flags().is_secondary()
    }

    /// 是否为QC失败; 对应flag: 0x200
    pub fn is_qc_fail(&self) -> bool {
        self.inner.flags().is_qc_fail()
    }

    /// 是否为重复reads; 对应flag: 0x400
    pub fn is_duplicate(&self) -> bool {
        self.inner.flags().is_duplicate()
    }

    /// 是否为补充比对信息; 对应flag: 0x800
    pub fn is_supplementary(&self) -> bool {
        self.inner.flags().is_supplementary()
    }
  
    /// 参考序列ID
    pub fn tid(&self) -> i32 {
        match self.inner.reference_sequence_id() {
            Some(Ok(id)) => id as i32,
            _ => -1,
        }
    }

    /// mate的参考序列ID
    pub fn mtid(&self) -> i32 {
        match self.inner.mate_reference_sequence_id() {
            Some(Ok(id)) => id as i32,
            _ => -1,
        }
    }


    pub fn insert_size(&self) -> i64 {
        self.inner.template_length() as i64
    }
}