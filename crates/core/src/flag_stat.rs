//! samtools flagstat 的实现

use bamqc_io::bam::{BamRecord};
use std::fmt;

#[derive(Debug)]
pub struct FlagStat {
    total: u64,
    primary: u64,
    secondary: u64,
    supplementary: u64,
    duplicate: u64,
    mapped: u64,
    primary_mapped: u64,
}



impl FlagStat {

    pub fn new() -> Self {

        Self {
            total: 0,
            primary: 0,
            secondary: 0,
            supplementary: 0,
            duplicate: 0,
            mapped: 0,
            primary_mapped: 0,
        }
    }
    pub fn update(&mut self, record: &BamRecord) {
        
        // QC失败的记录直接跳过，不计入任何统计
        if record.is_qc_fail() {
            return;
        }

        self.total += 1;

        // 按记录类型分类统计
        if !(record.is_secondary() | record.is_supplementary()) {
            self.primary += 1;
        }

    
        if record.is_secondary() {
            self.secondary += 1;
        }

        if record.is_supplementary() {
            self.supplementary += 1;
        }

        if record.is_duplicate() {
            self.duplicate += 1;
        }

        if !record.is_unmapped() {
            self.mapped += 1;
        }

        if !record.is_unmapped() & !(record.is_secondary() | record.is_supplementary()) {
            self.primary_mapped += 1;
        }


    }

    pub fn mapped_rate(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            self.mapped as f64 / self.total as f64
        }
    }
}

impl fmt::Display for FlagStat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // 简洁的统计信息输出格式
        writeln!(f, "total: {}", self.total)?;
        writeln!(f, "primary: {}", self.primary)?;
        writeln!(f, "secondary: {}", self.secondary)?;
        writeln!(f, "supplementary: {}", self.supplementary)?;
        writeln!(f, "duplicate: {}", self.duplicate)?;
        write!(f, "mapped: {} ({:.2}%)", self.mapped, self.mapped_rate() * 100.0)
    }
}