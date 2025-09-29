use clap::{Parser, Subcommand};
use bamqc_core::{
    PairOrientation, Strategy, compute_insert_size
};
use std::path::Path;
use std::fs::write;
use tracing::error;

/// BAM/CRAM文件质量控制工具组
#[derive(Parser)]
#[command(author, version, about = "BAM/CRAM文件质量控制工具组", long_about = None)]
struct Cli {
    /// 启用详细日志
    #[arg(short, long, global = true)]
    verbose: bool,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// 计算插入片段长度（与Picard CollectInsertSizeMetrics一致）
    InsertSize {
        /// 输入BAM/CRAM文件路径
        #[arg(short, long)]
        input: String,

        /// 输出文件路径（可选，如果不指定则输出到标准输出）
        #[arg(short, long)]
        output: Option<String>,

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
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    // 初始化日志
    let log_level = if cli.verbose { "debug" } else { "info" };
    tracing_subscriber::fmt()
        .with_env_filter(format!("bamqc={}", log_level))
        .init();

    match cli.command {
        Commands::InsertSize {
            input,
            output,
            include_duplicates,
            require_proper_pair,
            min_pct,
            pair_orientation,
            strategy,
        } => {
            handle_insert_size_command(
                &input,
                output,
                include_duplicates,
                require_proper_pair,
                min_pct,
                pair_orientation,
                strategy,
            )
        }
    }
}

/// 处理insert_size子命令
fn handle_insert_size_command(
    input: &str,
    output: Option<String>,
    include_duplicates: bool,
    require_proper_pair: bool,
    min_pct: f64,
    pair_orientation: PairOrientation,
    strategy: Strategy,
) -> Result<(), Box<dyn std::error::Error>> {
    // 验证输入文件存在
    if !Path::new(input).exists() {
        error!("输入文件不存在: {}", input);
        std::process::exit(1);
    }

    match compute_insert_size(
        input,
        include_duplicates,
        require_proper_pair,
        min_pct,
        pair_orientation,
        strategy,
    ) {
        Ok(median_size) => {
            let result = median_size.to_string();
            
            // 根据是否提供输出文件决定输出方式
            match output {
                Some(output_path) => {
                    // 写入文件
                    match write(&output_path, &result) {
                        Ok(_) => {
                            println!("结果已保存到文件: {}", output_path);
                            Ok(())
                        }
                        Err(e) => {
                            error!("写入文件失败 {}: {}", output_path, e);
                            std::process::exit(1);
                        }
                    }
                }
                None => {
                    // 输出到标准输出
                    println!("{}", result);
                    Ok(())
                }
            }
        }
        Err(e) => {
            error!("{}", e);
            std::process::exit(1);
        }
    }
}
