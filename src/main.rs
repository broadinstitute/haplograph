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
mod methyl;
mod eval;
mod extract;

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
    Haplograph {
        /// Input BAM file
        #[arg(short, long)]
        alignment_bam: String,
   
        /// Input FASTA file
        #[arg(short, long)]
        reference_fa: String,
        
        /// Sample ID
        #[arg(short, long)]
        sampleid: String,

        /// Output prefix
        #[arg(short, long, default_value = "./haplograph_")]
        output_prefix: String,
       
        /// a continuous genomic region as String (chromo:start-end) or a bed file as String
        #[arg(short, long)]
        locus: String,

        /// minimal variant allele frequency
        #[arg(short, long, default_value_t = 0.01)]
        var_frequency_min: f64,
        
        /// Minimal Supported Reads
        #[arg(short, long, default_value_t = 2)]
        min_reads: u8,

        ///window size
        #[arg(short, long, default_value_t = 100)]
        window_size: usize,
   
        ///if only primary reads are used
        #[arg(short, long, default_value = "false")]
        primary_only: bool, 

        ///output file format, accepted fasta, gfa, vcf
        #[arg(short, long, default_value = "gfa")]
        file_format: String, 

        /// Haplotype number
        #[arg(short, long, default_value_t = 2)]
        number_of_haplotypes: usize,

        /// heterozygous coverage fold threshold, > 3.0 is not heterozygous,  the smaller the more strict
        #[arg(short, long, default_value_t = 3.0)]
        coverage_fold_threshold: f64,

        /// methylation likelihood threshold, default to 0.5
        #[arg(short, long, default_value_t = 0.5)]
        threshold_methyl_likelihood: f32,

        /// Sequencing technology, accepted hifi, nanopore
        #[arg(short, long, default_value = "hifi")]
        detection_technology: String,

        /// Verbose output
        #[arg(long)]
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

        /// Germline only, default to false
        #[arg(short, long, default_value = "false")]
        major_haplotype_only: bool,

        /// Haplotype number
        #[arg(short, long, default_value_t = 2)]
        number_of_haplotypes: usize,

        /// heterozygous coverage fold threshold, > 3.0 is not heterozygous, the smaller the more strict
        #[arg(short, long, default_value_t = 3.0)]
        fold_threshold: f64,

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
        output_prefix: String,

        /// Sample ID
        #[arg(short, long)]
        sampleid: String,

        /// Reference FASTA file
        #[arg(short, long)]
        reference_fa: String,

        /// Maximum Haplotype number
        #[arg(short, long, default_value_t = 2)]
        maximum_haplotypes: usize,

        /// heterozygous coverage fold threshold, > 3.0 is not heterozygous,  the smaller the more strict
        #[arg(short, long, default_value_t = 3.0)]
        fold_threshold: f64,

        /// Sequencing technology, accepted hifi, nanopore
        #[arg(short, long, default_value = "hifi")]
        detection_technology: String,

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
    /// extract all the seqeunces from the bam file
    #[clap(arg_required_else_help = true)]
    Extract {
        /// Input Bam file
        #[arg(short, long)]
        bamfile: String,

        /// locus as String (chromo:start-end)
        #[arg(short, long)]
        locus: String,

        /// Output prefix
        #[arg(short, long, default_value = "haplograph_extract")]
        output_prefix: String,

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
        Commands::Haplograph {
            // Input BAM file
            alignment_bam,
            // Input FASTA file
            reference_fa,
            // Sample ID
            sampleid,
            // Output directory
            output_prefix,
            // either locus as String (chromo:start-end) or a bed file
            locus,
            // Limited size of the region
            var_frequency_min,
            // Minimal Supported Reads
            min_reads,
            //window size
            window_size,
            //if only primary reads are used
            primary_only, 
            //output file format
            file_format, 
            // haplotype number
            number_of_haplotypes,
            // heterozygous coverage fold threshold, > 3.0 is not heterozygous,  the smaller the more strict
            coverage_fold_threshold,
            // methylation likelihood threshold, default to 0.5
            threshold_methyl_likelihood,
            // Sequencing technology, accepted hifi, nanopore
            detection_technology,
            // Verbose output
            verbose,
        } => {
                // Validate format
                if file_format != "fasta" && file_format != "gfa" && file_format != "vcf" {
                    anyhow::bail!("Format must be either 'fasta' or 'gfa' or 'vcf', got: {}", file_format);
                }

                // Initialize logging
                if verbose {
                    std::env::set_var("RUST_LOG", "debug");
                } else {
                    std::env::set_var("RUST_LOG", "info");
                }
                env_logger::init();

                if locus.ends_with("bed.gz") || locus.ends_with(".bed") {
                    
                    info!("Starting Haplograph analysis");
                    info!("Input BAM: {}", alignment_bam);
                    info!("Region in {}", locus);
                    info!("Minimal vaf : {}", var_frequency_min);
                    info!("Minimal supported reads: {}", min_reads);
                    info!("Maximal Window size: {}", window_size);
                    info!("Primary only: {}", primary_only);
                    info!("Output prefix: {}", output_prefix);
                    info!("Default file format: {}", file_format);
                    info!("Verbose: {}", verbose);

                    let mut bam = util::open_bam_file(&alignment_bam);
                    let reference_seqs = util::get_all_ref_seq(&reference_fa);
                    // // Extract read sequences from BAM file using utility function
                    let bed_list = util::import_bed(&locus);
                    let mut windows = Vec::new();
                    for (chromosome, start, end) in bed_list {
                        for i in (start..end).step_by(window_size as usize) {
                            let end_pos = std::cmp::min(i + window_size , end);
                            windows.push((chromosome.clone(), i, end_pos));
                        }
                    }
                    windows.sort_by(|a, b| a.1.cmp(&b.1).then(a.2.cmp(&b.2)));
                    if file_format == "gfa" {
                        graph::start( &alignment_bam, &windows, &reference_seqs, &sampleid, min_reads as usize, threshold_methyl_likelihood, var_frequency_min, primary_only, &output_prefix, number_of_haplotypes)?;
                    } else {
                        hap::start( &alignment_bam.clone(), &windows, &reference_seqs, &sampleid, min_reads as usize, var_frequency_min, primary_only, &output_prefix, &file_format)?;
                    }

                } else {
                    let (chromosome, start, end) = util::split_locus(locus.clone());
                    info!("Starting Haplograph analysis");
                    info!("Input BAM: {}", alignment_bam);
                    info!("Region: {}:{}-{}", chromosome, start, end);
                    info!("Locus size: {}", end - start);
                    info!("Minimal vaf : {}", var_frequency_min);
                    info!("Minimal supported reads: {}", min_reads);
                    info!("Maximal Window size: {}", window_size);
                    info!("Primary read only: {}", primary_only);

                    let mut bam = util::open_bam_file(&alignment_bam);
                    let (reference_seqs, reference_chromosome_seqs) = util::get_ref_seq_from_chromosome(&reference_fa, &chromosome);
                    // // Extract read sequences from BAM file using utility function
                    let mut windows = Vec::new();
                    for i in (start..end).step_by(window_size as usize) {
                        let end_pos = std::cmp::min(i + window_size , end);
                        windows.push((chromosome.clone(), i, end_pos));
                    }
                    if file_format == "gfa" {
                        graph::start( &alignment_bam, &windows, &reference_chromosome_seqs, &sampleid, min_reads as usize, threshold_methyl_likelihood, var_frequency_min, primary_only, &output_prefix, number_of_haplotypes)?;
            
                    } else {
                        hap::start(&alignment_bam, &windows, &reference_chromosome_seqs, &sampleid, min_reads as usize, var_frequency_min, primary_only, &output_prefix, &file_format)?;
                    } 

                    let output_p = PathBuf::from(&output_prefix);
                    let graph_gfa = output_p.with_extension("gfa");
                    asm::start(&graph_gfa,  true, number_of_haplotypes,  &output_p, coverage_fold_threshold)?;                    
                    call::start(&graph_gfa, &reference_seqs, &sampleid, &output_prefix, number_of_haplotypes, coverage_fold_threshold, &detection_technology)?;

                }        

            }
        Commands::Assemble {
            graph_gfa,
            output_prefix,
            major_haplotype_only,
            number_of_haplotypes,
            fold_threshold,
            verbose,
        } => {
            // Initialize logging
            if verbose {
                std::env::set_var("RUST_LOG", "debug");
            } else {
                std::env::set_var("RUST_LOG", "info");
            }
            env_logger::init();
            asm::start(&graph_gfa,  major_haplotype_only, number_of_haplotypes,  &output_prefix, fold_threshold)?;
        }
        Commands::Call {
            gfa_file,
            output_prefix,
            sampleid,
            reference_fa,
            verbose,
            maximum_haplotypes,
            fold_threshold,
            detection_technology,
        } => {
            // Initialize logging
            if verbose {
                std::env::set_var("RUST_LOG", "debug");
            } else {
                std::env::set_var("RUST_LOG", "info");
            }
            env_logger::init();

            let reference_seqs = util::get_all_ref_seq(&reference_fa);
            call::start(&gfa_file, &reference_seqs, &sampleid, &output_prefix, maximum_haplotypes, fold_threshold, &detection_technology)?;
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
        Commands::Extract {
            bamfile,
            locus,
            output_prefix,
            verbose,
        } => {
            // Initialize logging
            if verbose {
                std::env::set_var("RUST_LOG", "debug");
            } else {
    }
                std::env::set_var("RUST_LOG", "info");
            env_logger::init();
            let (chromosome, start, end) = util::split_locus(locus.clone());
            let mut bam = util::open_bam_file(&bamfile);
            extract::start(&mut bam, &chromosome, start, end, false, output_prefix.clone().to_string())?;
        }
    }
    Ok(())
}
