use clap::Parser;
use log::{info};
use anyhow::{Result};
use std::path::PathBuf;

mod util;
mod hap;
mod intervals;
mod graph;

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
    output_prefix: String,
    
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

    ///if only primary reads are used
    #[arg(short, long, default_value = "false")]
    primary_only: bool, 

    ///output file format
    #[arg(short, long, default_value = "gfa")]
    default_file_format: String, 


    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

fn main() -> Result<()> {
    // Parse command line arguments
    let cli = Cli::parse();
    
    // Validate format
    if cli.default_file_format != "fasta" && cli.default_file_format != "gfa" {
        anyhow::bail!("Format must be either 'fasta' or 'gfa', got: {}", cli.default_file_format);
    }
    
    
    // Initialize logging
    if cli.verbose {
        std::env::set_var("RUST_LOG", "debug");
    } else {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::init();
    
    let (chromosome, start, end) = util::split_locus(cli.locus.clone());
    info!("Starting Haplograph analysis");
    info!("Input BAM: {}", cli.alignment_bam);
    info!("Region: {}:{}-{}", chromosome, start, end);
    info!("Locus size: {}", end - start);
    info!("Minimal vaf : {}", cli.frequency_min);
    info!("Minimal supported reads: {}", cli.min_reads);
    info!("Maximal Window size: {}", cli.window_size);
    info!("Primary only: {}", cli.primary_only);

    let mut bam = util::open_bam_file(&cli.alignment_bam);
    // // Extract read sequences from BAM file using utility function

    if cli.window_size  > end  - start  {
        let final_hap = intervals::start(&mut bam, &cli.reference_fa, &chromosome, start, end, &cli.sampleid, cli.min_reads as usize, cli.frequency_min, cli.primary_only, true)?;
        

    } else {
        let mut windows = Vec::new();
        for i in (start..end).step_by(cli.window_size as usize) {
            let end_pos = std::cmp::min(i + cli.window_size , end);
            windows.push((i, end_pos));
        }

        if cli.default_file_format == "fasta" {
            hap::start(&cli, &mut bam, &chromosome, &windows)?;

        } else if cli.default_file_format == "gfa" {
            graph::start( &cli, &mut bam, &windows)?;

        }
    }

    

    // let gfa_output = PathBuf::from(format!("{}/{}_{}_{}_{}_haplograph.gfa", cli.output, cli.sampleid, cli.chromosome, cli.start, cli.end));
    // let _ = build::start(&gfa_output, cli.kmer_size, reads, reference);
    // let reference_length = (cli.end - cli.start) as i64;
    // let haplotype_output = PathBuf::from(format!("{}/{}_{}_{}_{}_haplotypes.fasta", cli.output, cli.sampleid, cli.chromosome, cli.start, cli.end));
    // let _ = haplotypes::start(&gfa_output, reference_length, cli.bin_size, cli.pad_size, cli.min_read_ratio, cli.count_support, &haplotype_output, &cli.sampleid);
    
    Ok(())
}
