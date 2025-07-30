use clap::Parser;
use log::{info, error, warn};
use anyhow::{Result, Context};
use rust_htslib::bam::IndexedReader;
use std::path::{Path, PathBuf};
use url::Url;
use bio::io::fastq;
use rust_htslib::faidx::Reader;

mod util;
mod build;
mod haplotypes;
mod agg;

#[derive(Parser)]
#[command(name = "haplograph")]
#[command(about = "A bioinformatics tool for haplotype analysis")]
#[command(version)]
struct Cli {
    /// Input BAM file
    #[arg(short, long)]
    bam: String,

    /// Input FASTA file
    #[arg(short, long)]
    reference_fa: String,
    
    /// Sample ID
    #[arg(short, long)]
    sampleid: String,

    /// Output directory
    #[arg(short, long, default_value = ".")]
    output: String,

    /// kmer size to build graph
    #[arg(short, long, default_value = "21")]
    kmer_size: usize,
    
    /// Chromosome name
    #[arg(short, long)]
    chromosome: String,
    
    /// Start position
    #[arg(short, long)]
    start: u64,
    
    /// End position
    #[arg(short, long)]
    end: u64,
    
    /// Minimal Supported Reads
    #[arg(short, long)]
    min_reads: Option<u8>,

    /// Minimal mapping quality
    #[arg(short, long, default_value = "10")]
    quality: u8,
    
    /// bin size for haplotype reconstruction
    #[arg(short, long, default_value = "1000")]
    bin_size: i32, 

    ///pad size for adjacent window
    #[arg(short, long, default_value = "100")]
    pad_size:i64, 
    
    /// minimal read ratio, default 0.01
    #[arg(short, long, default_value = "0.01")]
    min_read_ratio:f64, 

    /// minimal count support default 2
    #[arg(short, long, default_value = "2")]
    count_support:usize, 

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

fn main() -> Result<()> {
    // Parse command line arguments
    let cli = Cli::parse();
    
    // Initialize logging
    if cli.verbose {
        std::env::set_var("RUST_LOG", "debug");
    } else {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::init();
    
    info!("Starting Haplograph analysis");
    info!("Input BAM: {}", cli.bam);
    info!("Region: {}:{}-{}", cli.chromosome, cli.start, cli.end);
    
    // Validate inputs
    let reads = if cli.start >= cli.end {
        anyhow::bail!("Start position must be less than end position");
    } else {
        // Continue with processing
        info!("Input validation passed");
        // Process the BAM file
        process_bam_file(&cli)
    };
    
    let reference = if cli.start >= cli.end {
        anyhow::bail!("Start position must be less than end position");
    } else {
        // Continue with processing
        info!("reference validation passed");
        // Process the BAM file
        process_fasta_file(&cli)
    };
    let gfa_output = PathBuf::from(format!("{}/{}_{}_{}_{}_haplograph.gfa", cli.output, cli.sampleid, cli.chromosome, cli.start, cli.end));
    let _ = build::start(&gfa_output, cli.kmer_size, reads, reference);
    let reference_length = (cli.end - cli.start) as i64;
    let haplotype_output = PathBuf::from(format!("{}/{}_{}_{}_{}_haplotypes.fasta", cli.output, cli.sampleid, cli.chromosome, cli.start, cli.end));
    let _ = haplotypes::start(&gfa_output, reference_length, cli.bin_size, cli.pad_size, cli.min_read_ratio, cli.count_support, &haplotype_output, &cli.sampleid);
    
    Ok(())
}

fn process_bam_file(cli: &Cli) -> Vec<fastq::Record> {
    
    // Open BAM file
    if cli.bam.to_string().starts_with("gs://") && std::env::var("GCS_OAUTH_TOKEN").is_err() {
        util::gcs_authorize_data_access();
    }
    let mut bam = if cli.bam.starts_with("gs://") {
        let url = Url::parse(&cli.bam).unwrap();
        IndexedReader::from_url(&url).unwrap()
    } else {
        IndexedReader::from_path(&cli.bam).unwrap()
    };
    
    info!("Successfully opened BAM file");
    
    // Extract reads from the specified region
    let reads = util::extract_reads_at_single_locus(
        &mut bam,
        &cli.chromosome,
        cli.start,
        cli.end,
        Some(cli.quality), // min_mapping_quality
    ).unwrap();
    
    info!("Extracted {} reads from region", reads.len());

    return reads;
}

fn process_fasta_file(cli: &Cli) -> Vec<fastq::Record> {
    
    // Open BAM file
    if cli.reference_fa.to_string().starts_with("gs://") && std::env::var("GCS_OAUTH_TOKEN").is_err() {
        util::gcs_authorize_data_access();
    }
    let mut fasta = if cli.reference_fa.starts_with("gs://") {
        let url = Url::parse(&cli.reference_fa).unwrap();
        Reader::from_url(&url).unwrap()
    } else {
        Reader::from_path(&cli.reference_fa).unwrap()
    };
    
    info!("Successfully opened FASTA file");
    
    // Extract reads from the specified region
    let reference_seqs = util::extract_fasta_seqs(
        &cli.sampleid,
        &mut fasta,
        &cli.chromosome,
        &cli.start,
        &cli.end,
        &"".to_string()
    ).unwrap();
    
    info!("Extracted {} sequences from region", reference_seqs.len());

   return reference_seqs;
}


fn write_fastq_output(reads: &[bio::io::fastq::Record], output_path: &str) -> Result<()> {
    use std::fs::File;
    use std::io::BufWriter;
    use bio::io::fastq::Writer;
    
    let file = File::create(output_path)
        .with_context(|| format!("Failed to create output file: {}", output_path))?;
    let writer = BufWriter::new(file);
    let mut fastq_writer = Writer::new(writer);
    
    for read in reads {
        fastq_writer.write_record(read)
            .with_context(|| "Failed to write FASTQ record")?;
    }
    
    Ok(())
}