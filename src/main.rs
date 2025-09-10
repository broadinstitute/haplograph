use log::{info};
use anyhow::{Result};
use clap::{Parser, Subcommand};
use std::path::PathBuf;

mod util;
mod hap;
mod intervals;
mod graph;
mod asm;
mod call;
mod eval;

#[derive(Parser)]
#[command(name = "haplograph")]
#[command(about = "A bioinformatics tool for haplotype analysis")]
#[command(version)]
pub struct Cli {
    #[clap(subcommand)]
    command: Commands,
 
}


#[derive(Debug, Subcommand)]

enum Commands {
    /// Extract Haplotypes from BAM file
    #[clap(arg_required_else_help = true)]
    Haplotype {
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
    },
    /// Assemble haplotypes from GFA file
    #[clap(arg_required_else_help = true)]
    Assemble {
        /// Input GFA file
        #[arg(short, long)]
        graph_gfa: PathBuf,
        
        /// Output prefix
        #[arg(short, long, default_value = "haplograph_asm")]
        output_prefix: PathBuf,

        /// Verbose output
       #[arg(short, long)]
       verbose: bool,
    },

    /// Call variants from VCF file
    #[clap(arg_required_else_help = true)]
    Call {
        /// Input VCF file
        #[arg(short, long)]
        gfa_file: PathBuf,

        /// Output prefix
        #[arg(short, long, default_value = "haplograph_call")]
        output_prefix: PathBuf,

        /// Sample ID
        #[arg(short, long)]
        sampleid: String,

        /// Reference FASTA file
        #[arg(short, long)]
        reference_fa: String,

        /// Verbose output
       #[arg(short, long)]
       verbose: bool,
    },
    /// Evaluate the accuracy of the haplotype calling
    #[clap(arg_required_else_help = true)]
    Evaluate {
        /// Input FASTA file
        #[arg(short, long)]
        truth_fasta: PathBuf,

        /// Input FASTA file
        #[arg(short, long)]
        query_fasta: PathBuf,

        /// Haplotype number
        #[arg(short, long, default_value = "2")]
        seq_number: usize,

        /// Output prefix
        #[arg(short, long, default_value = "haplograph_eval")]
        output_prefix: PathBuf,

        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },
 
}

fn main() -> Result<()> {
    // Parse command line arguments
    let cli = Cli::parse();

    let args = Cli::parse();
    match args.command {
        Commands::Haplotype {
            // Input BAM file
            alignment_bam,
            // Input FASTA file
            reference_fa,
            // Sample ID
            sampleid,
            // Output directory
            output_prefix,
            // Chromosome name
            locus,
            // Limited size of the region
            frequency_min,
            // Minimal Supported Reads
            min_reads,
            //window size
            window_size,
            //if only primary reads are used
            primary_only, 
            //output file format
            default_file_format, 
            // Verbose output
            verbose,
        } => {
                // Validate format
                if default_file_format != "fasta" && default_file_format != "gfa" {
                    anyhow::bail!("Format must be either 'fasta' or 'gfa', got: {}", default_file_format);
                }

                // Initialize logging
                if verbose {
                    std::env::set_var("RUST_LOG", "debug");
                } else {
                    std::env::set_var("RUST_LOG", "info");
                }
                env_logger::init();

                let (chromosome, start, end) = util::split_locus(locus.clone());
                info!("Starting Haplograph analysis");
                info!("Input BAM: {}", alignment_bam);
                info!("Region: {}:{}-{}", chromosome, start, end);
                info!("Locus size: {}", end - start);
                info!("Minimal vaf : {}", frequency_min);
                info!("Minimal supported reads: {}", min_reads);
                info!("Maximal Window size: {}", window_size);
                info!("Primary only: {}", primary_only);

                let mut bam = util::open_bam_file(&alignment_bam);
                let reference_seqs = util::get_chromosome_ref_seq(&reference_fa, &chromosome);
                // // Extract read sequences from BAM file using utility function
                let mut windows = Vec::new();
                for i in (start..end).step_by(window_size as usize) {
                    let end_pos = std::cmp::min(i + window_size , end);
                    windows.push((i, end_pos));
                }
        
                if default_file_format == "fasta" {
                    hap::start(&mut bam, &chromosome, &windows, locus, &reference_seqs, &sampleid, min_reads as usize, frequency_min, primary_only, &output_prefix)?;
        
                } else if default_file_format == "gfa" {
                    graph::start( &mut bam, &windows, &locus, &reference_seqs, &sampleid, min_reads as usize, frequency_min, primary_only, &output_prefix)?;
        
                }
            }
        Commands::Assemble {
            graph_gfa,
            output_prefix,
            verbose,
        } => {
            // Initialize logging
            if verbose {
                std::env::set_var("RUST_LOG", "debug");
            } else {
                std::env::set_var("RUST_LOG", "info");
            }
            env_logger::init();
            asm::start(&graph_gfa, &output_prefix)?;
        }
        Commands::Call {
            gfa_file,
            output_prefix,
            sampleid,
            reference_fa,
            verbose,
        } => {
            // Initialize logging
            if verbose {
                std::env::set_var("RUST_LOG", "debug");
            } else {
                std::env::set_var("RUST_LOG", "info");
            }
            env_logger::init();

            let reference_seqs = util::get_all_ref_seq(&reference_fa);
            call::start(&gfa_file, &reference_seqs, &sampleid, &output_prefix)?;
        }
        Commands::Evaluate {
            truth_fasta,
            query_fasta,
            seq_number,
            output_prefix,
            verbose,
        } => {

            // Initialize logging
            if verbose {
                std::env::set_var("RUST_LOG", "debug");
            } else {
                std::env::set_var("RUST_LOG", "info");
            }
            env_logger::init();
            eval::start(&truth_fasta, &query_fasta, seq_number, &output_prefix)?;
        }
    }
     
    Ok(())
}
