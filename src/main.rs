use anyhow::Result;
use bio::io::fasta::Reader as FastaReader;
use clap::{Args, Parser, Subcommand};
use log::{info, warn};
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

mod asm;
mod call;
mod eval;
mod extract;
mod graph;
mod hap;
mod haplopan;
mod intervals;
mod merge;
mod methyl;
mod util;

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

    /// Merge per-window haplograph outputs into full-length FASTA, VCF, and methylation BED
    #[command(name = "merge", arg_required_else_help = true)]
    MergeOutputs {
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

        /// Overlap between adjacent sub-loci in bp (must match the haplograph run)
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
    reference_fa: String,
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
    overlap_bp: usize,
    /// When set, a chunk whose window graph is too complex is split into smaller
    /// overlapping sub-chunks, assembled independently, and merged back together.
    adaptive_subdivide: bool,
}

// Adaptive subdivision thresholds/limits (only used when `adaptive_subdivide`).
const ADAPTIVE_MAX_DEPTH: usize = 3;
const ADAPTIVE_MIN_SIZE: usize = 50_000;
const ADAPTIVE_NODE_THRESHOLD: usize = 20_000;
const ADAPTIVE_PARALLEL_THRESHOLD: usize = 32;

/// Complexity proxy: too many nodes, or an interval with too many parallel nodes
/// (SV bubbles / high diversity), which blows up het identification and the DP.
fn graph_is_complex(node_info: &HashMap<String, asm::NodeInfo>) -> bool {
    if node_info.len() > ADAPTIVE_NODE_THRESHOLD {
        return true;
    }
    let parallel = asm::find_parallele_nodes(node_info);
    let max_parallel = parallel.values().map(|s| s.len()).max().unwrap_or(0);
    max_parallel > ADAPTIVE_PARALLEL_THRESHOLD
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
    run_haplograph_segment_depth(config, chromosome, start, end, segment_prefix, 0)
}

fn run_haplograph_segment_depth(
    config: &HaplographRunConfig,
    chromosome: &str,
    start: usize,
    end: usize,
    segment_prefix: &str,
    depth: usize,
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
    // graph::start writes `{prefix}.gfa` by appending; use the same construction
    // (with_extension would corrupt prefixes that contain a '.', e.g. `hla.r0_0`).
    let graph_gfa = PathBuf::from(format!("{segment_prefix}.gfa"));

    // Adaptive divide-and-conquer: only when enabled, and only if this graph is
    // genuinely complex, split into smaller overlapping sub-chunks and merge.
    if config.adaptive_subdivide
        && depth < ADAPTIVE_MAX_DEPTH
        && end.saturating_sub(start) > 2 * ADAPTIVE_MIN_SIZE
    {
        let (node_info, _edge_info) = asm::load_graph(&graph_gfa)?;
        if graph_is_complex(&node_info) {
            let sub_size = (end.saturating_sub(start) / 2).max(ADAPTIVE_MIN_SIZE);
            let sub_chunks =
                merge::plan_locus_chunks(chromosome, start, end, sub_size, config.overlap_bp)?;
            if sub_chunks.len() > 1 {
                info!(
                    "Adaptive subdivide {chromosome}:{start}-{end} (depth {depth}, {} nodes) into {} sub-chunks",
                    node_info.len(),
                    sub_chunks.len()
                );
                for (i, (c, s, e)) in sub_chunks.iter().enumerate() {
                    run_haplograph_segment_depth(
                        config,
                        c,
                        *s,
                        *e,
                        &format!("{segment_prefix}_{i}"),
                        depth + 1,
                    )?;
                }
                let locus_str = format!("{chromosome}:{start}-{end}");
                merge::start(
                    segment_prefix,
                    &locus_str,
                    config.number_of_haplotypes,
                    &config.sampleid,
                    &config.reference_fa,
                    sub_size,
                    config.overlap_bp,
                )?;
                return Ok(());
            }
        }
    }

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

/// Choose a per-region chunk size. With `adaptive` enabled, large regions are
/// split into more, smaller chunks (better parallel load balancing) down to a
/// floor, while never exceeding the user-requested size. `window_size` is
/// deliberately NOT scaled here — it is a resolution/correctness knob, not a
/// parallelization knob (see the roadmap's window_size caveat).
fn adaptive_locus_size(region_len: usize, requested: usize, adaptive: bool) -> usize {
    const MIN_CHUNK: usize = 50_000;
    const TARGET_CHUNKS: usize = 256;
    if !adaptive || region_len <= requested || requested == 0 {
        return requested;
    }
    let by_target = region_len / TARGET_CHUNKS;
    by_target.clamp(MIN_CHUNK, requested)
}

/// Genome-wide gather: index each per-region VCF and `bcftools concat` them into
/// a single sorted output. Failures (e.g. bcftools missing) are logged, not
/// fatal — the per-region VCFs remain usable.
fn concat_region_vcfs(kind: &str, region_prefixes: &[String], output_prefix: &str) -> Result<()> {
    use std::process::Command;

    let inputs: Vec<String> = region_prefixes
        .iter()
        .map(|p| format!("{p}.{kind}.vcf.gz"))
        .filter(|f| std::path::Path::new(f).exists())
        .collect();
    if inputs.is_empty() {
        warn!("No {kind} VCFs found to concatenate");
        return Ok(());
    }

    for f in &inputs {
        let ok = Command::new("bcftools")
            .args(["index", "-f", "-t", f])
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if !ok {
            warn!("bcftools index failed for {f}; skipping genome-wide {kind} concat");
            return Ok(());
        }
    }

    let out = format!("{output_prefix}.{kind}.vcf.gz");
    let mut cmd = Command::new("bcftools");
    cmd.args(["concat", "-a", "-Oz", "-o", &out]);
    cmd.args(&inputs);
    match cmd.status() {
        Ok(s) if s.success() => {
            let _ = Command::new("bcftools").args(["index", "-t", &out]).status();
            info!(
                "Wrote genome-wide {kind} VCF: {out} ({} regions)",
                inputs.len()
            );
        }
        _ => warn!("bcftools concat failed for {kind}; per-region {kind} VCFs are still available"),
    }
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

        /// Recursively subdivide a sub-locus when its window graph is too complex
        #[arg(long, default_value_t = false)]
        adaptive_subdivide: bool,

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

        /// Genomic locus to analyze (chromo:start-end)
        #[arg(short, long)]
        locus: String,

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

        /// Number of haplotypes to select from the pangenome (1 = haploid, 2 = diploid, etc.)
        #[arg(long, default_value_t = 2)]
        ploidy: usize,

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

    /// Whole-genome / multi-region analysis from a BED file: plans overlapping chunks per
    /// region, assembles/calls all chunks in parallel, stitches each region, and concatenates
    /// genome-wide germline/somatic VCFs.
    #[clap(arg_required_else_help = true)]
    Wgs {
        /// Input BAM file
        #[arg(short, long)]
        alignment_bam: String,

        /// Input Reference FASTA file
        #[arg(short, long)]
        reference_fa: String,

        /// Sample ID
        #[arg(short, long)]
        sampleid: String,

        /// Output prefix (per-region outputs are `{prefix}.rN.*`, genome-wide are `{prefix}.*`)
        #[arg(short, long, default_value = "./haplograph_wgs")]
        output_prefix: String,

        /// BED file of regions to analyze (whole genome = one line per contig)
        #[arg(short, long)]
        bed_file: String,

        /// minimal variant allele frequency
        #[arg(short, long, default_value_t = 0.01)]
        var_frequency_min: f64,

        /// Minimal Supported Reads
        #[arg(short, long, default_value_t = 2)]
        min_reads: u8,

        /// window size (resolution knob; not scaled by region length)
        #[arg(short, long, default_value_t = 100)]
        window_size: usize,

        /// if only primary reads are used
        #[arg(short, long, default_value_t = false)]
        primary_only: bool,

        /// if use pileup to extract reads
        #[arg(long, default_value_t = false)]
        pileup: bool,

        /// Haplotype number, currently only support 1 or 2
        #[arg(short, long, default_value_t = 2)]
        number_of_haplotypes: usize,

        /// methylation likelihood threshold, default to 0.5
        #[arg(short, long, default_value_t = 0.5)]
        threshold_methyl_likelihood: f32,

        /// Sequencing technology, accepted hifi, nanopore
        #[arg(short, long, default_value = "hifi")]
        detection_technology: String,

        /// maximal sub-locus (chunk) size
        #[arg(long, default_value_t = 200_000)]
        maximal_locus_size: usize,

        /// Overlap between adjacent chunks (bp)
        #[arg(long, default_value_t = merge::DEFAULT_OVERLAP_BP)]
        overlap_bp: usize,

        /// Scale chunk size down for very large regions (more, smaller chunks for load balancing)
        #[arg(long, default_value_t = false)]
        adaptive_chunking: bool,

        /// Skip chunks whose mean read depth is below this threshold
        #[arg(long, default_value_t = 2.0)]
        min_mean_depth: f64,

        /// Cap how many chunks assemble concurrently to bound peak memory (0 = unbounded).
        /// Inner per-window parallelism still fills idle cores within each chunk.
        #[arg(long, default_value_t = 0)]
        max_concurrent_chunks: usize,

        /// Recursively subdivide a chunk when its window graph is too complex
        #[arg(long, default_value_t = false)]
        adaptive_subdivide: bool,

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
            adaptive_subdivide,
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
                reference_fa: reference_fa.clone(),
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
                overlap_bp,
                adaptive_subdivide,
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

                let processed: Result<()> =
                    eligible_loci
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
                } else {
                    info!("Merging haplotypes");
                    merge::start(
                        &output_prefix,
                        &locus,
                        number_of_haplotypes,
                        &sampleid,
                        &reference_fa,
                        maximal_locus_size,
                        overlap_bp,
                    )?;
                    info!("Merging haplotypes completed");
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
                    run_haplograph_segment(&run_config, &chromosome, start, end, &output_prefix)?;
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
            locus,
            output_prefix,
            rollingkmer_list,
            sample_id,
            data_technology,
            pileup,
            //window size
            window_size,
            ploidy,
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

            let rollingkmer_list_vec: Vec<usize> = rollingkmer_list
                .split(",")
                .map(|x| x.parse::<usize>().unwrap())
                .collect();
            let (selected_haplotypes, _distance) = haplopan::start(
                &alignment_bam,
                &pangenome_fasta,
                &locus,
                &rollingkmer_list_vec,
                &output_prefix.to_string(),
                &data_technology.to_string(),
                &sample_id,
                ploidy,
            )?;
            let is_homozygous = {
                let unique_count = selected_haplotypes
                    .iter()
                    .collect::<std::collections::HashSet<_>>()
                    .len();
                unique_count < selected_haplotypes.len()
            };

            let tmp_bam_path = PathBuf::from(format!("{}.tmp.sorted.bam", output_prefix));
            let tmp_fasta_path = PathBuf::from(format!("{}.tmp.ref.fasta", output_prefix));

            let reference_seqs = util::get_all_ref_seq(&tmp_fasta_path.display().to_string());



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
                    asm::start(&graph_gfa, number_of_haplotypes, &output_prefix)?;
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
                DevToolsCommands::MergeOutputs {
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
            eval::start(
                &truth_fasta,
                &query_fasta,
                seq_number,
                &output_prefix,
                as_genotyper,
            )?;
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
            let mut bam = util::open_bam_file(&bamfile)?;
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
        Commands::Wgs {
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
            number_of_haplotypes,
            threshold_methyl_likelihood,
            detection_technology,
            maximal_locus_size,
            overlap_bp,
            adaptive_chunking,
            min_mean_depth,
            max_concurrent_chunks,
            adaptive_subdivide,
            verbose,
        } => {
            env_logger::Builder::from_default_env()
                .filter_level(if verbose {
                    log::LevelFilter::Debug
                } else {
                    log::LevelFilter::Info
                })
                .init();

            info!("Starting whole-genome Haplograph analysis");
            info!("Input BAM: {}", alignment_bam);
            info!("Regions from BED: {}", bed_file);
            info!("Output prefix: {}", output_prefix);
            info!("Chunk size: {} (adaptive={})", maximal_locus_size, adaptive_chunking);
            info!("Chunk overlap: {} bp", overlap_bp);

            if maximal_locus_size <= overlap_bp {
                anyhow::bail!(
                    "maximal_locus_size ({}) must be greater than overlap ({} bp)",
                    maximal_locus_size,
                    overlap_bp
                );
            }

            let bed_list = util::import_bed(&bed_file);
            if bed_list.is_empty() {
                anyhow::bail!("No regions found in BED file {}", bed_file);
            }

            if util::is_remote_bam(&alignment_bam) {
                warn!(
                    "Alignment BAM '{}' looks remote (gs://, s3://, http(s)://). Random seeks under \
                     high chunk-level parallelism are slow and error-prone; consider localizing the \
                     BAM (or a per-contig slice) before running whole-genome analysis.",
                    alignment_bam
                );
            }

            // Load reference once and index a single record per contig, so each
            // per-chromosome config carries only that contig (not the whole genome).
            let all_records = util::get_all_ref_seq(&reference_fa);
            let mut chrom_record: HashMap<String, bio::io::fastq::Record> = HashMap::new();
            for rec in &all_records {
                chrom_record
                    .entry(rec.id().to_string())
                    .or_insert_with(|| rec.clone());
            }
            drop(all_records);

            let mut configs: HashMap<String, Arc<HaplographRunConfig>> = HashMap::new();
            for (chrom, _, _) in &bed_list {
                if configs.contains_key(chrom) {
                    continue;
                }
                let rec = chrom_record.get(chrom).cloned().ok_or_else(|| {
                    anyhow::anyhow!("chromosome {chrom} not found in reference {reference_fa}")
                })?;
                let single = vec![rec];
                configs.insert(
                    chrom.clone(),
                    Arc::new(HaplographRunConfig {
                        alignment_bam: alignment_bam.clone(),
                        reference_fa: reference_fa.clone(),
                        reference_seqs: single.clone(),
                        reference_chromosome_seqs: single,
                        sampleid: sampleid.clone(),
                        min_reads,
                        window_size,
                        threshold_methyl_likelihood,
                        var_frequency_min,
                        primary_only,
                        pileup,
                        number_of_haplotypes,
                        detection_technology: detection_technology.clone(),
                        overlap_bp,
                        adaptive_subdivide,
                    }),
                );
            }

            // Plan regions and flatten all chunks into a single job list; the
            // chunk is the primary parallel axis (many chunks fill the cores).
            struct RegionInfo {
                chrom: String,
                start: usize,
                end: usize,
                prefix: String,
                size: usize,
            }
            let mut regions: Vec<RegionInfo> = Vec::new();
            let mut jobs: Vec<(Arc<HaplographRunConfig>, String, usize, usize, String)> = Vec::new();
            for (region_idx, (chrom, start, end)) in bed_list.iter().enumerate() {
                let region_len = end.saturating_sub(*start);
                if region_len == 0 {
                    continue;
                }
                let size = adaptive_locus_size(region_len, maximal_locus_size, adaptive_chunking);
                let region_prefix = format!("{output_prefix}.r{region_idx}");
                let chunks = merge::plan_locus_chunks(chrom, *start, *end, size, overlap_bp)?;
                let cfg = configs.get(chrom).unwrap().clone();
                for (ci, (c, s, e)) in chunks.iter().enumerate() {
                    let seg_prefix = format!("{region_prefix}_{ci}");
                    jobs.push((cfg.clone(), c.clone(), *s, *e, seg_prefix));
                }
                regions.push(RegionInfo {
                    chrom: chrom.clone(),
                    start: *start,
                    end: *end,
                    prefix: region_prefix,
                    size,
                });
            }
            info!(
                "Whole-genome plan: {} regions, {} chunks",
                regions.len(),
                jobs.len()
            );

            // Assemble + call every chunk. The chunk is the primary parallel axis;
            // a single chunk failure is logged and skipped rather than aborting the
            // whole genome. When `max_concurrent_chunks` is set, chunks are processed
            // in bounded batches to cap peak memory while inner per-window
            // parallelism still keeps idle cores busy within each chunk.
            let process_chunk = |(cfg, chrom, s, e, seg_prefix): &(
                Arc<HaplographRunConfig>,
                String,
                usize,
                usize,
                String,
            )| {
                let mean_depth = util::with_thread_local_bam(&cfg.alignment_bam, |bam| {
                    util::mean_depth(bam, chrom, *s, *e).ok()
                });
                if let Some(d) = mean_depth {
                    if d < min_mean_depth {
                        warn!(
                            "Skipping chunk {chrom}:{s}-{e} (mean depth {:.2} < {:.2})",
                            d, min_mean_depth
                        );
                        return Ok(());
                    }
                }
                // asm.rs can abort via panic (`.unwrap()`, `panic!`) on hard
                // graphs, not only via `Result::Err`. A panic inside a rayon
                // parallel iterator unwinds and aborts the whole parallel region,
                // so contain it here: one bad chunk is logged and skipped.
                let r = match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    run_haplograph_segment(cfg, chrom, *s, *e, seg_prefix)
                })) {
                    Ok(res) => res,
                    Err(payload) => {
                        let msg = payload
                            .downcast_ref::<&str>()
                            .map(|s| s.to_string())
                            .or_else(|| payload.downcast_ref::<String>().cloned())
                            .unwrap_or_else(|| "unknown panic".to_string());
                        Err(anyhow::anyhow!("panicked: {msg}"))
                    }
                };
                if let Err(err) = &r {
                    warn!("Chunk {chrom}:{s}-{e} failed: {err}");
                }
                r
            };

            let results: Vec<Result<()>> = if max_concurrent_chunks == 0 {
                jobs.par_iter().map(process_chunk).collect()
            } else {
                let mut acc = Vec::with_capacity(jobs.len());
                for batch in jobs.chunks(max_concurrent_chunks) {
                    acc.extend(batch.par_iter().map(process_chunk).collect::<Vec<_>>());
                }
                acc
            };
            let failed = results.iter().filter(|r| r.is_err()).count();
            if failed > 0 {
                warn!("{failed}/{} chunks failed (see logs above)", jobs.len());
            }

            // Stitch each region (in parallel across regions). A merge panic in
            // one region must not abort the whole genome, mirroring the
            // chunk-level containment above.
            regions.par_iter().for_each(|region| {
                let locus_str = format!("{}:{}-{}", region.chrom, region.start, region.end);
                let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    merge::start(
                        &region.prefix,
                        &locus_str,
                        number_of_haplotypes,
                        &sampleid,
                        &reference_fa,
                        region.size,
                        overlap_bp,
                    )
                }));
                match outcome {
                    Ok(Ok(())) => {}
                    Ok(Err(err)) => warn!("Merge failed for region {locus_str}: {err}"),
                    Err(_) => warn!("Merge panicked for region {locus_str}; skipping"),
                }
            });

            // Genome-wide gather in genomic order.
            let mut ordered = regions.iter().collect::<Vec<_>>();
            ordered.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
            let region_prefixes: Vec<String> =
                ordered.iter().map(|r| r.prefix.clone()).collect();
            concat_region_vcfs("germline", &region_prefixes, &output_prefix)?;
            concat_region_vcfs("somatic", &region_prefixes, &output_prefix)?;
            info!("Whole-genome Haplograph analysis complete");
        }
    }
    Ok(())
}
