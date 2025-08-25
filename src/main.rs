use clap::Parser;
use log::{info};
use anyhow::{Result};

mod util;
mod hap;

#[derive(Parser)]
#[command(name = "haplograph")]
#[command(about = "A bioinformatics tool for haplotype analysis")]
#[command(version)]
pub struct Cli {
    /// Input BAM file
    #[arg(short, long)]
    alignment_bam: String,

    /// Input FASTA file
    #[arg(short, long)]
    reference_fa: String,
    
    /// Sample ID
    #[arg(short, long)]
    sampleid: String,

    /// Output directory
    #[arg(short, long, default_value = ".")]
    output: String,
    
    /// Chromosome name
    #[arg(short, long)]
    locus: String,

    /// Limited size of the region
    #[arg(short, long, default_value = "0.01")]
    frequency_min: f64,
    
    /// Minimal Supported Reads
    #[arg(short, long, default_value = "2")]
    min_reads: u8,

    /// Minimal mapping quality
    #[arg(short, long, default_value = "10")]
    quality: u8,

    ///window size
    #[arg(short, long, default_value = "1000")]
    window_size: usize,

    ///pad size for adjacent window
    #[arg(short, long, default_value = "100")]
    pad_size:i64, 


    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

fn main() -> Result<()> {
    // Parse command line arguments
    let cli = Cli::parse();
    let (chromosome, start, end) = util::split_locus(cli.locus.clone());
    
    // Initialize logging
    if cli.verbose {
        std::env::set_var("RUST_LOG", "debug");
    } else {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::init();
    
    info!("Starting Haplograph analysis");
    info!("Input BAM: {}", cli.alignment_bam);
    info!("Region: {}:{}-{}", chromosome, start, end);
    info!("Locus size: {}", end - start);
    info!("Minimal vaf : {}", cli.frequency_min);
    info!("Minimal supported reads: {}", cli.min_reads);

    info!("Window size: {}", cli.window_size);
    info!("Pad size: {}", cli.pad_size);

    if cli.window_size  > end  - start  {
        // if the window size is larger than the region size, call the function directly
        hap::start(&cli, &chromosome, start as usize, end as usize)?;
    } else {
        let mut windows = Vec::new();
        let window_num = (end - start) / cli.window_size ;
        for i in (start..end).step_by(cli.window_size as usize) {
            let end_pos = std::cmp::min(i + cli.window_size , end);
            windows.push((i, end_pos));
        }
        for window in windows {
            hap::start(&cli, &chromosome, window.0 as usize, window.1 as usize)?;
        }
    }

    

    // let gfa_output = PathBuf::from(format!("{}/{}_{}_{}_{}_haplograph.gfa", cli.output, cli.sampleid, cli.chromosome, cli.start, cli.end));
    // let _ = build::start(&gfa_output, cli.kmer_size, reads, reference);
    // let reference_length = (cli.end - cli.start) as i64;
    // let haplotype_output = PathBuf::from(format!("{}/{}_{}_{}_{}_haplotypes.fasta", cli.output, cli.sampleid, cli.chromosome, cli.start, cli.end));
    // let _ = haplotypes::start(&gfa_output, reference_length, cli.bin_size, cli.pad_size, cli.min_read_ratio, cli.count_support, &haplotype_output, &cli.sampleid);
    
    Ok(())
}
