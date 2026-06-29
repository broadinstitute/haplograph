use anyhow::Result;
use clap::{Args, Parser, Subcommand};
use log::{info, warn};
use rayon::prelude::*;
use std::path::PathBuf;
use bio::io::fasta::Reader as FastaReader;
use std::collections::HashMap;
use std::sync::Arc;


mod asm;
mod call;
mod eval;
mod extract;
mod graph;
mod hap;
mod intervals;
mod merge;
mod methyl;
mod util;
mod haplopan;

#[derive(Debug, Subcommand)]
enum DevToolsCommands {
    /// Assemble haplotypes from GFA file
    #[clap(arg_required_else_help = true)]
    Assemble {
        /// Input GFA file
        #[arg(short, long)]
        graph_gfa: PathBuf,

        /// Output prefix
        #[arg(short, long, default_value = "haplograph_asm")]
        output_prefix: PathBuf,

        /// Haplotype number
        #[arg(short, long, default_value_t = 2)]
        number_of_haplotypes: usize,

        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },

    /// Call variants from GFA file
    #[clap(arg_required_else_help = true)]
    Call {
        /// Input GFA file
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

        /// Sequencing technology, accepted hifi, nanopore
        #[arg(short, long, default_value = "hifi")]
        detection_technology: String,

        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },
}



#[derive(Parser)]
#[command(name = "haplograph")]
#[command(about = "A bioinformatics tool for haplotype analysis")]
#[command(version)]
pub struct Cli {
    #[command(flatten)]
    global: GlobalOpts,

    #[clap(subcommand)]
    command: Commands,
}

#[derive(Args, Debug)]
pub struct GlobalOpts {
    /// Rayon worker threads for parallel sub-loci / windows (overrides RAYON_NUM_THREADS)
    #[arg(long, global = true)]
    pub threads: Option<usize>,
}

const MINIMAL_GAP_LENGTH: usize = 50;

struct HaplographRunConfig {
    alignment_bam: String,
    reference_seqs: Vec<bio::io::fastq::Record>,
    reference_chromosome_seqs: Vec<bio::io::fastq::Record>,
    sampleid: String,
    min_reads: u8,
    window_size: usize,
    threshold_methyl_likelihood: f32,
    var_frequency_min: f64,
    primary_only: bool,
    pileup: bool,
    number_of_haplotypes: usize,
    detection_technology: String,
}

fn build_windows(
    chromosome: &str,
    start: usize,
    end: usize,
    window_size: usize,
) -> Vec<(String, usize, usize)> {
    let mut windows = Vec::new();
    for i in (start..end).step_by(window_size) {
        let end_pos = std::cmp::min(i + window_size, end);
        windows.push((chromosome.to_string(), i, end_pos));
    }
    windows
}

fn run_haplograph_segment(
    config: &HaplographRunConfig,
    chromosome: &str,
    start: usize,
    end: usize,
    segment_prefix: &str,
) -> Result<()> {
    let windows = build_windows(chromosome, start, end, config.window_size);
    graph::start(
        &config.alignment_bam,
        &windows,
        &config.reference_chromosome_seqs,
        &config.sampleid,
        config.min_reads as usize,
        config.threshold_methyl_likelihood,
        config.var_frequency_min,
        config.primary_only,
        config.pileup,
        &segment_prefix.to_string(),
        MINIMAL_GAP_LENGTH,
    )?;

    let output_p = PathBuf::from(segment_prefix);
    let graph_gfa = output_p.with_extension("gfa");
    let (primary_haplotypes, node_info, edge_info) =
        asm::start(&graph_gfa, config.number_of_haplotypes, &output_p)?;
    call::start(
        &graph_gfa,
        &config.reference_seqs,
        &primary_haplotypes,
        &config.sampleid,
        &segment_prefix.to_string(),
        config.number_of_haplotypes,
        &config.detection_technology,
        Some((node_info, edge_info)),
    )?;
    Ok(())
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Haplograph Analysis from BAM file for a continuous genomic region (usually a locus > 1kb)
    #[clap(arg_required_else_help = true)]
    Haplograph {
        /// Input BAM file
        #[arg(short, long)]
        alignment_bam: String,

        /// Input Reference FASTA file
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
        #[arg(short, long, default_value_t = false)]
        primary_only: bool,

        ///if use pileup to extract reads
        #[arg(long, default_value_t = false)]
        pileup: bool,

        ///output file format, accepted fasta, gfa, vcf
        #[arg(short, long, default_value = "gfa")]
        file_format: String,

        /// Haplotype number, currently only support 1 or 2
        #[arg(short, long, default_value_t = 2)]
        number_of_haplotypes: usize,

        /// methylation likelihood threshold, default to 0.5
        #[arg(short, long, default_value_t = 0.5)]
        threshold_methyl_likelihood: f32,

        /// Sequencing technology, accepted hifi, nanopore
        #[arg(short, long, default_value = "hifi")]
        detection_technology: String,

        /// maxial locus size
        #[arg(long, default_value_t = 200000)]
        maximal_locus_size: usize,

        /// Overlap between adjacent sub-loci when splitting large regions (bp)
        #[arg(long, default_value_t = merge::DEFAULT_OVERLAP_BP)]
        overlap_bp: usize,

        /// Skip sub-loci whose mean read depth is below this threshold
        #[arg(long, default_value_t = 2.0)]
        min_mean_depth: f64,

        /// Verbose output
        #[arg(long)]
        verbose: bool,
    },
    /// Somatic-aware Haplotype Interval Analysis from BAM file for a set of genomic region defined by a bed file(usually a locus < 1kb, e.g. STRs)
    #[clap(arg_required_else_help = true)]
    Haplointervals {
        /// Input BAM file
        #[arg(short, long)]
        alignment_bam: String,

        /// Input ReferenceFASTA file
        #[arg(short, long)]
        reference_fa: String,

        /// Sample ID
        #[arg(short, long)]
        sampleid: String,

        /// Output prefix
        #[arg(short, long, default_value = "./haplograph_")]
        output_prefix: String,

        /// a bed file containing the genomic regions to be analyzed
        #[arg(short, long)]
        bed_file: String,

        /// minimal variant allele frequency
        #[arg(short, long, default_value_t = 0.01)]
        var_frequency_min: f64,

        /// Minimal Supported Reads
        #[arg(short, long, default_value_t = 2)]
        min_reads: u8,

        ///Maximalwindow size
        #[arg(short, long, default_value_t = 1000)]
        window_size: usize,

        ///if only primary reads are used
        #[arg(short, long, default_value_t = false)]
        primary_only: bool,

        ///if use pileup to extract reads
        #[arg(long, default_value_t = false)]
        pileup: bool,

        /// methylation likelihood threshold, default to 0.5
        #[arg(short, long, default_value_t = 0.5)]
        threshold_methyl_likelihood: f32,

        /// Verbose output
        #[arg(long)]
        verbose: bool,
    },

    /// HaploPan Analysis leveraging pangenomes to resolve complex haplotypes (e.g. duplications, deletions, translocations, etc.)
    #[clap(arg_required_else_help = true)]
    Haplopan {
        /// Input Pangenome FASTA file
        #[arg(short, long)]
        pangenome_fasta: PathBuf,

        /// Input BAM file
        #[arg(short, long)]
        alignment_bam: PathBuf,

        /// Output prefix
        #[arg(short, long, default_value = "haplograph_pangenome")]
        output_prefix: String,

        /// Rolling kmer list
        #[arg(short, long, default_value = "31")]
        rollingkmer_list: String,

        /// Sample ID
        #[arg(short, long)]
        sample_id: String,

        /// Sequencing technology, accepted hifi, nanopore, sr
        #[arg(short, long, default_value = "hifi")]
        data_technology: String,
        //window size
        #[arg(short, long, default_value_t = 100)]
        window_size: usize,

        ///if use pileup to extract reads
        #[arg(long, default_value_t = false)]
        pileup: bool,

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

        /// Repurpose as genotyper, default to false
        #[arg(short, long, default_value = "false")]
        as_genotyper: bool,

        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },
    /// extract all the seqeunces in an intervalfrom the bam file
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

        /// if use pileup to extract reads
        #[arg(short, long, default_value = "false")]
        pileup: bool,

        /// Sample ID
        #[arg(short, long)]
        sampleid: String,

        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },
    /// Development tools
    #[clap(subcommand)]
    DevTools(DevToolsCommands),

    /// Merge per-window haplograph outputs into full-length FASTA, VCF, and methylation BED
    #[clap(arg_required_else_help = true)]
    Merge {
        /// Base output prefix used for the large-locus haplograph run (segments are `{prefix}_0`, `{prefix}_1`, …)
        #[arg(short, long)]
        output_prefix: String,

        /// Full genomic locus (chrom:start-end) that was analyzed
        #[arg(short, long)]
        locus: String,

        /// Input Reference FASTA file
        #[arg(short, long)]
        reference_fa: String,

        /// Sample ID
        #[arg(short, long)]
        sampleid: String,

        /// Haplotype number, currently only support 1 or 2
        #[arg(short, long, default_value_t = 2)]
        number_of_haplotypes: usize,

        /// Overlap between adjacent sub-loci in bp (must match the haplograph run; default 500)
        #[arg(long, default_value_t = merge::DEFAULT_OVERLAP_BP)]
        overlap_bp: usize,

        /// Maximal sub-locus size used during haplograph (must match the haplograph run; default 200000)
        #[arg(long, default_value_t = 200_000)]
        maximal_locus_size: usize,

        /// Verbose output
        #[arg(long)]
        verbose: bool,
    },
}

fn main() -> Result<()> {
    let args = Cli::parse();
    util::init_rayon_threads(args.global.threads)?;
    match args.command {
        Commands::Haplograph {
            alignment_bam,
            reference_fa,
            sampleid,
            output_prefix,
            locus,
            var_frequency_min,
            min_reads,
            window_size,
            primary_only,
            pileup,
            file_format,
            number_of_haplotypes,
            threshold_methyl_likelihood,
            detection_technology,
            maximal_locus_size,
            overlap_bp,
            min_mean_depth,
            verbose,
        } => {
            // Validate format
            if file_format != "fasta" && file_format != "gfa" && file_format != "vcf" {
                anyhow::bail!(
                    "Format must be either 'fasta' or 'gfa' or 'vcf', got: {}",
                    file_format
                );
            }

            // Initialize logging
            env_logger::Builder::from_default_env()
                .filter_level(if verbose {
                    log::LevelFilter::Debug
                } else {
                    log::LevelFilter::Info
                })
                .init();

            let (chromosome, start, end) = util::split_locus(locus.clone());
            info!("Starting Haplograph analysis");
            info!("Input BAM: {}", alignment_bam);
            info!("Region: {}:{}-{}", chromosome, start, end);
            info!("Locus size: {}", end - start);
            info!("Minimal vaf : {}", var_frequency_min);
            info!("Minimal supported reads: {}", min_reads);
            info!("Maximal Window size: {}", window_size);
            info!("Primary read only: {}", primary_only);
            info!("Minimum mean depth for sub-loci: {}", min_mean_depth);
            info!("Sub-locus overlap: {} bp", overlap_bp);

            let (reference_seqs, reference_chromosome_seqs) =
                util::get_ref_seq_from_chromosome(&reference_fa, &chromosome);
            let run_config = Arc::new(HaplographRunConfig {
                alignment_bam: alignment_bam.clone(),
                reference_seqs,
                reference_chromosome_seqs,
                sampleid: sampleid.clone(),
                min_reads,
                window_size,
                threshold_methyl_likelihood,
                var_frequency_min,
                primary_only,
                pileup,
                number_of_haplotypes,
                detection_technology: detection_technology.clone(),
            });

            if end - start > maximal_locus_size {
                info!("Locus size is larger than the maximal locus size: {}, splitting into overlapping intervals", end - start);
                if maximal_locus_size <= overlap_bp {
                    anyhow::bail!(
                        "maximal_locus_size ({}) must be greater than overlap ({} bp)",
                        maximal_locus_size,
                        overlap_bp
                    );
                }
                let step = maximal_locus_size - overlap_bp;
                let mut locus_list = Vec::new();
                // Windows of length maximal_locus_size, advanced by step so adjacent windows share overlap_bp.
                let mut start_pos = start;
                while start_pos < end {
                    let end_pos = std::cmp::min(start_pos + maximal_locus_size, end);
                    locus_list.push((chromosome.clone(), start_pos, end_pos));
                    if end_pos >= end {
                        break;
                    }
                    start_pos += step;
                }
                locus_list.sort_by(|a, b| a.1.cmp(&b.1).then(a.2.cmp(&b.2)));

                let eligible_loci: Vec<(String, usize, usize)> = locus_list
                    .par_iter()
                    .filter_map(|locus| {
                        let mean_depth = util::with_thread_local_bam(&alignment_bam, |bam| {
                            util::mean_depth(bam, &locus.0, locus.1, locus.2).ok()
                        })?;
                        if mean_depth < min_mean_depth {
                            warn!(
                                "Skipping sub-locus {}:{}-{} (mean depth {:.2} < {:.2})",
                                locus.0, locus.1, locus.2, mean_depth, min_mean_depth
                            );
                            return None;
                        }
                        Some(locus.clone())
                    })
                    .collect();

                let processed: Result<()> = eligible_loci
                    .par_iter()
                    .enumerate()
                    .try_for_each(|(output_index, locus)| {
                        info!("Interval: {}:{}-{}", locus.0, locus.1, locus.2);
                        let segment_prefix = format!("{}_{}", output_prefix, output_index);
                        run_haplograph_segment(
                            &run_config,
                            &locus.0,
                            locus.1,
                            locus.2,
                            &segment_prefix,
                        )
                    });
                processed?;

                if eligible_loci.is_empty() {
                    warn!("All sub-loci were skipped due to low coverage");
                }
            } else {
                let mean_depth = util::with_thread_local_bam(&alignment_bam, |bam| {
                    util::mean_depth(bam, &chromosome, start, end)
                })?;
                if mean_depth < min_mean_depth {
                    warn!(
                        "Skipping locus {}:{}-{} (mean depth {:.2} < {:.2})",
                        chromosome, start, end, mean_depth, min_mean_depth
                    );
                } else {
                    run_haplograph_segment(
                        &run_config,
                        &chromosome,
                        start,
                        end,
                        &output_prefix,
                    )?;
                }
            }
        }
        Commands::Haplointervals {
            alignment_bam,
            reference_fa,
            sampleid,
            output_prefix,
            bed_file,
            var_frequency_min,
            min_reads,
            window_size,
            primary_only,
            pileup,
            threshold_methyl_likelihood,
            verbose,
        } => {
            // Initialize logging
            env_logger::Builder::from_default_env()
                .filter_level(if verbose {
                    log::LevelFilter::Debug
                } else {
                    log::LevelFilter::Info
                })
                .init();

            info!("Starting Haplograph analysis");
            info!("Input BAM: {}", alignment_bam);
            info!("Region in {}", bed_file);
            info!("Minimal vaf : {}", var_frequency_min);
            info!("Minimal supported reads: {}", min_reads);
            info!("Maximal Window size: {}", window_size);
            info!("Primary only: {}", primary_only);
            info!("Output prefix: {}", output_prefix);
            info!("Default output file format: {}", "vcf");
            info!("Verbose: {}", verbose);

            let reference_seqs = util::get_all_ref_seq(&reference_fa);
            let bed_list = util::import_bed(&bed_file);
            let mut windows = Vec::new();
            for (chromosome, start, end) in bed_list {
                if end - start > window_size {
                    info!(
                        "Window size is too large, skipping region: {}:{}-{}",
                        chromosome, start, end
                    );
                    continue;
                } else {
                    for i in (start..end).step_by(window_size) {
                        let end_pos = std::cmp::min(i + window_size, end);
                        windows.push((chromosome.clone(), i, end_pos));
                    }
                }
            }
            windows.sort_by(|a, b| a.1.cmp(&b.1).then(a.2.cmp(&b.2)));

            hap::start(
                &alignment_bam.clone(),
                &windows,
                &reference_seqs,
                &sampleid,
                min_reads as usize,
                var_frequency_min,
                primary_only,
                pileup,
                &output_prefix,
                &"vcf".to_string(),
                threshold_methyl_likelihood,
            )?;
        }

        Commands::Haplopan {
            pangenome_fasta,
            alignment_bam,
            output_prefix,
            rollingkmer_list,
            sample_id,
            data_technology,
            pileup,
            //window size
            window_size,

            verbose,
        } => {

            // Initialize logging
            env_logger::Builder::from_default_env()
            .filter_level(if verbose {
                log::LevelFilter::Debug
            } else {
                log::LevelFilter::Info
            })
            .init();
            info!("Starting HaploPan analysis");
            info!("Input BAM: {}", alignment_bam.display());
            info!("Input Pangenome FASTA: {}", pangenome_fasta.display());
            info!("Rolling kmer list: {}", rollingkmer_list);
            info!("Output prefix: {}", output_prefix);
            info!("Verbose: {}", verbose);

            let rollingkmer_list_vec = rollingkmer_list.split(",").map(|x| x.parse::<usize>().unwrap()).collect();
            haplopan::start(&alignment_bam, &pangenome_fasta, &rollingkmer_list_vec, &output_prefix.to_string(), &data_technology.to_string(), &sample_id)?;

            let tmp_bam_path = PathBuf::from(format!("{}.tmp.sorted.bam", output_prefix));
            let tmp_fasta_path = PathBuf::from(format!("{}.tmp.ref.fasta", output_prefix));

            let reference_seqs = util::get_all_ref_seq(&tmp_fasta_path.display().to_string());
            
            let mut final_fasta_seq = HashMap::new();
            for record in reference_seqs.iter(){
                let mut windows = Vec::new();
                let start = 0;
                let end = record.seq().len();
                let chromosome = record.id().to_string();
                for i in (start..end).step_by(window_size) {
                    let end_pos = std::cmp::min(i + window_size, end);
                    windows.push((chromosome.clone(), i, end_pos));                 
                }
                graph::start(
                    &tmp_bam_path.display().to_string(),
                    &windows,
                    &reference_seqs,
                    &sample_id,
                    1,
                    0.5,
                    0.0,
                    false,
                    pileup,
                    &format!("{}_{}_tmp", output_prefix, chromosome),
                    MINIMAL_GAP_LENGTH
                )?;
                let output_p = PathBuf::from(&format!("{}_{}_tmp", output_prefix, chromosome));
                let graph_gfa = output_p.with_extension("gfa");
                asm::start(
                    &graph_gfa,
                    1,
                    &output_p,
                )?;
                let fasta_reader = FastaReader::from_file(format!("{}.fasta", &output_p.display().to_string()))?;
                for record in fasta_reader.records() {
                    let record = record.expect("Failed to read FASTA record");
                    final_fasta_seq.insert(record.id().to_string(), String::from_utf8_lossy(record.seq()).to_string());
                }
            }
            util::write_fasta(&final_fasta_seq, &PathBuf::from(format!("{}.final.fasta", output_prefix)))?;
        }

        Commands::DevTools(dev_tools_cmd) => {
            match dev_tools_cmd {
                DevToolsCommands::Assemble {
                    graph_gfa,
                    output_prefix,
                    number_of_haplotypes,
                    verbose,
                } => {
                    // Initialize logging
                    env_logger::Builder::from_default_env()
                        .filter_level(if verbose {
                            log::LevelFilter::Debug
                        } else {
                            log::LevelFilter::Info
                        })
                        .init();
                    asm::start(
                        &graph_gfa,
                        number_of_haplotypes,
                        &output_prefix,
                    )?;
                }
                DevToolsCommands::Call {
                    gfa_file,
                    output_prefix,
                    sampleid,
                    reference_fa,
                    verbose,
                    maximum_haplotypes,
                    detection_technology,
                } => {
                    // Initialize logging
                    env_logger::Builder::from_default_env()
                        .filter_level(if verbose {
                            log::LevelFilter::Debug
                        } else {
                            log::LevelFilter::Info
                        })
                        .init();

                    let reference_seqs = util::get_all_ref_seq(&reference_fa);
                    let (primary_haplotypes, node_info, edge_info) = asm::start(
                        &gfa_file,
                        maximum_haplotypes,
                        &PathBuf::from(&output_prefix.clone()),
                    )?;
                    call::start(
                        &gfa_file,
                        &reference_seqs,
                        &primary_haplotypes,
                        &sampleid,
                        &output_prefix,
                        maximum_haplotypes,
                        &detection_technology,
                        Some((node_info, edge_info)),
                    )?;
                    
                }
            }
        }
        Commands::Evaluate {
            truth_fasta,
            query_fasta,
            seq_number,
            output_prefix,
            as_genotyper,
            verbose,
        } => {
            // Initialize logging
            env_logger::Builder::from_default_env()
                .filter_level(if verbose {
                    log::LevelFilter::Debug
                } else {
                    log::LevelFilter::Info
                })
                .init();
            eval::start(&truth_fasta, &query_fasta, seq_number, &output_prefix, as_genotyper)?;
        }
        Commands::Extract {
            bamfile,
            locus,
            output_prefix,
            pileup,
            sampleid,
            verbose,
        } => {
            // Initialize logging
            env_logger::Builder::from_default_env()
                .filter_level(if verbose {
                    log::LevelFilter::Debug
                } else {
                    log::LevelFilter::Info
                })
                .init();
            let (chromosome, start, end) = util::split_locus(locus.clone());
            let mut bam = util::open_bam_file(&bamfile);
            extract::start(
                &mut bam,
                &chromosome,
                start,
                end,
                false,
                output_prefix.clone().to_string(),
                sampleid.clone(),
                pileup,
            )?;
        }
        Commands::Merge {
            output_prefix,
            locus,
            reference_fa,
            sampleid,
            number_of_haplotypes,
            overlap_bp,
            maximal_locus_size,
            verbose,
        } => {
            env_logger::Builder::from_default_env()
                .filter_level(if verbose {
                    log::LevelFilter::Debug
                } else {
                    log::LevelFilter::Info
                })
                .init();
            merge::start(
                &output_prefix,
                &locus,
                number_of_haplotypes,
                &sampleid,
                &reference_fa,
                maximal_locus_size,
                overlap_bp,
            )?;
        }
    }
    Ok(())
}
