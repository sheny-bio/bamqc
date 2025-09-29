use clap::Parser;
use bamqc_core::{
    PairOrientation, Strategy, compute_insert_size
};
use std::path::Path;
use tracing::error;

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
