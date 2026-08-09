use anyhow::{bail, Context, Result as AnyhowResult};
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation;
use bio::io::fasta::Reader as FastaReader;
use bio::io::fastq;
use flate2::read::GzDecoder;
use indicatif::ProgressBar;
use log::debug;
use log::info;
use log::warn;
use ndarray::s;
use ndarray::{Array1, Array2};
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rayon::prelude::*;
use rust_htslib::bam::{self, record::Aux, IndexedReader, Read as BamRead};
use rust_htslib::faidx;
use statrs::distribution::ContinuousCDF;
use statrs::distribution::Normal;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::path::PathBuf;
use url::Url;

pub fn reverse_complement(kmer: &str) -> String {
    kmer.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            _ => panic!("Unexpected character: {}", c),
        })
        .collect()
}

pub fn write_fasta(
    all_sequences: &HashMap<String, String>,
    output_filename: &PathBuf,
) -> AnyhowResult<()> {
    if let Some(parent) = output_filename.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut file = File::create(output_filename)?;
    let chars_per_line = 60;
    for (header, sequence) in all_sequences.iter() {
        writeln!(file, ">{}", header)?;
        let seq_len = sequence.len();
        let full_lines = seq_len / chars_per_line;
        for i in 0..full_lines {
            let start = i * chars_per_line;
            let end = start + chars_per_line;
            writeln!(file, "{}", &sequence[start..end])?;
        }
        if seq_len % chars_per_line != 0 {
            writeln!(file, "{}", &sequence[full_lines * chars_per_line..])?;
        }
    }
    Ok(())
}

pub fn gcs_gcloud_is_installed() -> bool {
    // Check if gcloud is installed on the PATH
    // Suppress stdout and stderr to prevent them from printing to the screen
    let mut cmd = std::process::Command::new("gcloud");
    cmd.arg("version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .is_ok()
}

pub fn gcs_authorize_data_access() -> AnyhowResult<()> {
    // Check if gcloud is installed on the PATH
    if !gcs_gcloud_is_installed() {
        bail!("gcloud is not installed on the PATH (needed to authorize gs:// access)");
    }

    // Execute the command and capture the output
    let output = std::process::Command::new("gcloud")
        .args(["auth", "application-default", "print-access-token"])
        .output()
        .context("failed to execute `gcloud auth application-default print-access-token`")?;

    if !output.status.success() {
        bail!(
            "gcloud auth failed: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        );
    }

    // Decode the output and remove trailing newline
    let token = String::from_utf8(output.stdout)
        .context("gcloud returned a non-UTF8 access token")?
        .trim_end()
        .to_string();

    // Set the environment variable
    unsafe {
        std::env::set_var("GCS_OAUTH_TOKEN", token);
    }
    Ok(())
}

// Function to get a mapping between read group and sample name from a BAM header.
pub fn get_rg_to_sm_mapping(bam: &IndexedReader) -> HashMap<String, String> {
    let header = bam::Header::from_template(bam.header());

    let rg_sm_map: HashMap<String, String> = header
        .to_hashmap()
        .into_iter()
        .flat_map(|(_, records)| records)
        .filter(|record| record.contains_key("ID") && record.contains_key("SM"))
        .map(|record| (record["ID"].clone(), record["SM"].clone()))
        .collect();

    rg_sm_map
}

pub fn get_sm_name_from_rg(
    read: &bam::Record,
    rg_sm_map: &HashMap<String, String>,
) -> AnyhowResult<String> {
    let rg = read.aux(b"RG")?;

    if let Aux::String(v) = rg {
        if let Some(sm) = rg_sm_map.get(v) {
            Ok(sm.to_owned())
        } else {
            Err(anyhow::anyhow!(
                "Sample name not found for read group: {}",
                v
            ))
        }
    } else {
        Err(anyhow::anyhow!("Read group is not a string"))
    }
}

/// Step size for sampled mean-depth estimation on large intervals.
pub const MEAN_DEPTH_SAMPLE_STEP: usize = 1024;

/// Configure the global rayon thread pool (`--threads`). Must run before any `par_iter`.
/// When `threads` is `None`, rayon uses `RAYON_NUM_THREADS` or the machine default.
pub fn init_rayon_threads(threads: Option<usize>) -> AnyhowResult<()> {
    let Some(n) = threads else {
        return Ok(());
    };
    anyhow::ensure!(n >= 1, "--threads must be at least 1 (got {n})");
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build_global()
        .context("failed to initialize rayon thread pool (already initialized?)")?;
    info!("Rayon thread pool: {n} worker(s)");
    Ok(())
}

/// A BAM path is "remote" (random seeks are expensive under high concurrency)
/// if it is a cloud/URL path rather than a local file.
pub fn is_remote_bam(bam_path: &str) -> bool {
    let lower = bam_path.to_ascii_lowercase();
    lower.starts_with("gs://")
        || lower.starts_with("s3://")
        || lower.starts_with("http://")
        || lower.starts_with("https://")
        || lower.starts_with("ftp://")
}

thread_local! {
    static THREAD_BAM: RefCell<Option<(String, IndexedReader)>> = const { RefCell::new(None) };
}

/// Reuse one `IndexedReader` per rayon worker thread instead of opening per window.
pub fn with_thread_local_bam<F, R>(bam_path: &str, f: F) -> R
where
    F: FnOnce(&mut IndexedReader) -> R,
{
    THREAD_BAM.with(|cell| {
        let mut slot = cell.borrow_mut();
        let need_open = slot
            .as_ref()
            .map(|(path, _)| path != bam_path)
            .unwrap_or(true);
        if need_open {
            // A failed open is fatal for every caller (they all need the BAM), so
            // surface the full actionable error chain instead of an opaque
            // `unwrap()` panic. Transient/auth failures are already retried inside
            // `open_bam_file`, so reaching here means the BAM is genuinely unusable.
            let reader = open_bam_file(bam_path).unwrap_or_else(|err| panic!("{err:#}"));
            *slot = Some((bam_path.to_string(), reader));
        }
        f(&mut slot.as_mut().unwrap().1)
    })
}

/// How many times to attempt opening a *remote* BAM before giving up. Remote
/// opens fail transiently under high concurrency and when a GCS token expires
/// mid-run; retrying (with a re-auth) recovers both without aborting the run.
const REMOTE_BAM_OPEN_ATTEMPTS: usize = 3;

pub fn open_bam_file(alignment_bam: &str) -> AnyhowResult<IndexedReader> {
    if is_remote_bam(alignment_bam) {
        open_remote_bam(alignment_bam)
    } else {
        let reader = IndexedReader::from_path(alignment_bam).with_context(|| {
            format!(
                "Failed to open local BAM '{alignment_bam}' \
                 (is it indexed? expected '{alignment_bam}.bai')"
            )
        })?;
        debug!("Successfully opened BAM file");
        Ok(reader)
    }
}

/// Point `CURL_CA_BUNDLE` at a CA bundle that actually exists on this system, so
/// htslib/libcurl can do TLS. Only set it when unset AND a real bundle is found:
/// pointing the variable at a missing path (e.g. the Debian bundle on macOS)
/// breaks not just libcurl but every child process that reads it — including the
/// `gcloud` we shell out to for GCS auth.
fn ensure_curl_ca_bundle() {
    if std::env::var_os("CURL_CA_BUNDLE").is_some() {
        return;
    }
    // Ordered by prevalence; the first existing path wins.
    const CANDIDATES: &[&str] = &[
        "/etc/ssl/certs/ca-certificates.crt",   // Debian / Ubuntu
        "/etc/pki/tls/certs/ca-bundle.crt",     // RHEL / CentOS / Fedora
        "/etc/ssl/ca-bundle.pem",               // openSUSE
        "/etc/ssl/cert.pem",                    // macOS / Alpine
        "/opt/homebrew/etc/openssl@3/cert.pem", // macOS Homebrew (Apple silicon)
        "/usr/local/etc/openssl@3/cert.pem",    // macOS Homebrew (Intel)
    ];
    match CANDIDATES.iter().find(|p| Path::new(p).exists()) {
        Some(path) => {
            unsafe {
                std::env::set_var("CURL_CA_BUNDLE", path);
            }
            debug!("Set CURL_CA_BUNDLE to {path}");
        }
        None => {
            warn!(
                "No CA certificate bundle found in standard locations; leaving CURL_CA_BUNDLE \
                 unset. If remote BAM opens fail with a TLS error, set CURL_CA_BUNDLE (and/or \
                 SSL_CERT_FILE) to your system's CA bundle."
            );
        }
    }
}

fn open_remote_bam(alignment_bam: &str) -> AnyhowResult<IndexedReader> {
    // GCS needs an OAuth token; fetch one up front if we don't already have it.
    if alignment_bam.starts_with("gs://") && std::env::var("GCS_OAUTH_TOKEN").is_err() {
        gcs_authorize_data_access()
            .context("failed to obtain a GCS access token before opening the remote BAM")?;
    }
    // htslib/libcurl needs a CA bundle for TLS
    // (https://github.com/rust-bio/rust-htslib/issues/404#issuecomment-1905264507).
    ensure_curl_ca_bundle();
    let url =
        Url::parse(alignment_bam).with_context(|| format!("Invalid BAM URL '{alignment_bam}'"))?;

    let mut last_err = None;
    for attempt in 1..=REMOTE_BAM_OPEN_ATTEMPTS {
        match IndexedReader::from_url(&url) {
            Ok(reader) => {
                debug!("Successfully opened remote BAM file (attempt {attempt})");
                return Ok(reader);
            }
            Err(err) => {
                warn!(
                    "Opening remote BAM '{alignment_bam}' failed \
                     (attempt {attempt}/{REMOTE_BAM_OPEN_ATTEMPTS}): {err:?}"
                );
                // A common cause under long/parallel runs is an expired GCS token;
                // refresh it before the next attempt.
                if alignment_bam.starts_with("gs://") {
                    if let Err(auth_err) = gcs_authorize_data_access() {
                        warn!("GCS re-authorization failed: {auth_err:#}");
                    }
                }
                last_err = Some(err);
                if attempt < REMOTE_BAM_OPEN_ATTEMPTS {
                    std::thread::sleep(std::time::Duration::from_millis(300 * attempt as u64));
                }
            }
        }
    }

    Err(last_err.expect("loop ran at least once")).with_context(|| {
        format!(
            "Failed to open remote BAM '{alignment_bam}' after {REMOTE_BAM_OPEN_ATTEMPTS} attempts. \
             Check that (1) the index '{alignment_bam}.bai' exists and is readable, \
             (2) GCS credentials are valid (`gcloud auth application-default login`), and \
             (3) htslib was built with libcurl/GCS support."
        )
    })
}

/// Mean read depth over `[start, end)` on `chromosome`, counting uncovered positions as zero.
/// Large intervals use position sampling every [`MEAN_DEPTH_SAMPLE_STEP`] bp.
pub fn mean_depth(
    bam: &mut IndexedReader,
    chromosome: &str,
    start: usize,
    end: usize,
) -> AnyhowResult<f64> {
    let region_len = end.saturating_sub(start);
    if region_len == 0 {
        return Ok(0.0);
    }
    if bam.header().tid(chromosome.as_bytes()).is_none() {
        anyhow::bail!("Chromosome {chromosome} not found in BAM header");
    }
    bam.fetch((chromosome.as_bytes(), start as u64, end as u64))?;

    let mean = if region_len <= MEAN_DEPTH_SAMPLE_STEP * 2 {
        let mut depth_at = vec![0u32; region_len];
        for pileup_result in bam.pileup() {
            let pileup = pileup_result?;
            let pos = pileup.pos() as usize;
            if pos >= start && pos < end {
                depth_at[pos - start] = pileup.depth();
            }
        }
        depth_at.iter().map(|d| *d as u64).sum::<u64>() as f64 / region_len as f64
    } else {
        let sample_count = region_len.div_ceil(MEAN_DEPTH_SAMPLE_STEP);
        let mut depth_samples = vec![0u32; sample_count];
        for pileup_result in bam.pileup() {
            let pileup = pileup_result?;
            let pos = pileup.pos() as usize;
            if pos < start || pos >= end {
                continue;
            }
            let offset = pos - start;
            if offset.is_multiple_of(MEAN_DEPTH_SAMPLE_STEP) {
                depth_samples[offset / MEAN_DEPTH_SAMPLE_STEP] = pileup.depth();
            }
        }
        depth_samples.iter().map(|d| *d as u64).sum::<u64>() as f64 / sample_count as f64
    };

    info!(
        "Mean depth for {}:{}-{} = {:.2}x",
        chromosome, start, end, mean
    );
    Ok(mean)
}

fn read_fasta_streaming(reference_fa: &str) -> Vec<fastq::Record> {
    debug!("Opening FASTA file (streaming): {}", reference_fa);
    let file = File::open(reference_fa).expect("Failed to open FASTA file");
    let reader: Box<dyn BufRead> = if reference_fa.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let fasta_reader = FastaReader::new(reader);
    let mut reference_seqs = Vec::new();
    for result in fasta_reader.records() {
        let record = result.expect("Failed to read FASTA record");
        let seq_id = record.id().to_string();
        let sequence = String::from_utf8_lossy(record.seq()).to_string();
        reference_seqs.push(fastq::Record::with_attrs(
            &seq_id,
            None,
            sequence.as_bytes(),
            vec![30; sequence.len()].as_slice(),
        ));
    }
    reference_seqs
}

/// FASTA accessor backed by faidx (one chromosome at a time).
pub struct ReferenceFasta {
    reader: faidx::Reader,
}

impl ReferenceFasta {
    pub fn open(reference_fa: &str) -> AnyhowResult<Self> {
        if reference_fa.ends_with(".gz") {
            anyhow::bail!(
                "indexed FASTA access requires an uncompressed FASTA with a .fai index; got {reference_fa}"
            );
        }
        let path = Path::new(reference_fa);
        let fai_path = format!("{reference_fa}.fai");
        if !Path::new(&fai_path).exists() {
            faidx::build(path)
                .map_err(|err| anyhow::anyhow!("failed to build FASTA index: {err}"))?;
        }
        let reader = faidx::Reader::from_path(path)
            .map_err(|err| anyhow::anyhow!("failed to open indexed FASTA {reference_fa}: {err}"))?;
        Ok(Self { reader })
    }

    pub fn sequence_names(&self) -> AnyhowResult<Vec<String>> {
        self.reader
            .seq_names()
            .map_err(|err| anyhow::anyhow!("failed to read FASTA sequence names: {err}"))
    }

    pub fn fetch_chromosome_record(&self, chromosome: &str) -> AnyhowResult<fastq::Record> {
        let seq_len = self.reader.fetch_seq_len(chromosome) as usize;
        let sequence = self
            .reader
            .fetch_seq_string(chromosome, 0, seq_len)
            .map_err(|err| anyhow::anyhow!("failed to fetch {chromosome}: {err}"))?;
        Ok(fastq::Record::with_attrs(
            chromosome,
            None,
            sequence.as_bytes(),
            vec![30; sequence.len()].as_slice(),
        ))
    }

    pub fn fetch_all_chromosome_records(&self) -> AnyhowResult<Vec<fastq::Record>> {
        self.sequence_names()?
            .into_iter()
            .map(|chrom| self.fetch_chromosome_record(&chrom))
            .collect()
    }
}

pub fn get_all_ref_seq(reference_fa: &String) -> Vec<fastq::Record> {
    if reference_fa.ends_with(".gz") {
        return read_fasta_streaming(reference_fa);
    }
    match ReferenceFasta::open(reference_fa) {
        Ok(genome) => genome
            .fetch_all_chromosome_records()
            .unwrap_or_else(|_| read_fasta_streaming(reference_fa)),
        Err(_) => read_fasta_streaming(reference_fa),
    }
}

pub fn get_ref_seq_from_chromosome(
    reference_fa: &String,
    chromosome: &str,
) -> (Vec<fastq::Record>, Vec<fastq::Record>) {
    if reference_fa.ends_with(".gz") {
        let reference_seqs = get_all_ref_seq(reference_fa);
        let reference_seqs_chromosome = reference_seqs
            .iter()
            .filter(|r| r.id() == chromosome)
            .cloned()
            .collect();
        return (reference_seqs, reference_seqs_chromosome);
    }

    let genome = ReferenceFasta::open(reference_fa)
        .unwrap_or_else(|_| panic!("Failed to open reference FASTA: {reference_fa}"));
    let reference_seqs = genome
        .fetch_all_chromosome_records()
        .unwrap_or_else(|_| get_all_ref_seq(reference_fa));
    let reference_seqs_chromosome = vec![genome
        .fetch_chromosome_record(chromosome)
        .unwrap_or_else(|_| {
            reference_seqs
                .iter()
                .find(|r| r.id() == chromosome)
                .cloned()
                .expect("chromosome not found in reference")
        })];
    (reference_seqs, reference_seqs_chromosome)
}

pub fn split_locus(locus: String) -> (String, usize, usize) {
    let parts: Vec<&str> = locus.split(':').collect();
    let chromosome = parts[0].to_string();
    let start: usize = parts[1]
        .split("-")
        .collect::<Vec<&str>>()
        .first()
        .unwrap()
        .parse()
        .unwrap();
    let end: usize = parts[1]
        .split("-")
        .collect::<Vec<&str>>()
        .last()
        .unwrap()
        .parse()
        .unwrap();
    (chromosome, start, end)
}

/// Parses a run-length CIGAR string ("23=4I1X20=") into `(length, op)` pairs.
pub fn parse_cigar(cigar: &str) -> Vec<(usize, char)> {
    let mut ops = Vec::new();
    let mut num = String::new();
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num.push(c);
        } else {
            ops.push((num.parse().expect("cigar length"), c));
            num.clear();
        }
    }
    ops
}

/// Rewrites a CIGAR to eliminate a length-1 `X` sitting directly against an
/// `I` or `D` run whenever the X's "orphaned" base (the one with no partner
/// on the gapped side) also occurs inside that run. In that case the run can
/// be losslessly re-anchored so the X becomes a plain match and the indel is
/// redistributed around it: same total ref/alt content, no substitution.
/// Left untouched when the orphaned base doesn't occur in the adjacent run -
/// that's very likely a genuine substitution, not a homopolymer-registration
/// artifact.
pub fn normalize_homopolymer_indels(cigar: &str, ref_seq: &str, alt_seq: &str) -> String {
    let ops = parse_cigar(cigar);
    let ref_bytes = ref_seq.as_bytes();
    let alt_bytes = alt_seq.as_bytes();

    let mut result: Vec<(usize, char)> = Vec::new();
    let mut ref_pos = 0usize;
    let mut alt_pos = 0usize;
    let mut i = 0usize;

    let push_op = |result: &mut Vec<(usize, char)>, len: usize, op: char| {
        if len == 0 {
            return;
        }
        if let Some(last) = result.last_mut() {
            if last.1 == op {
                last.0 += len;
                return;
            }
        }
        result.push((len, op));
    };

    while i < ops.len() {
        let (len, op) = ops[i];

        // I(n) immediately followed by a single mismatch: the mismatch's ref
        // base may really belong inside the insertion run.
        if op == 'I' && i + 1 < ops.len() && ops[i + 1] == (1, 'X') {
            let ins = &alt_bytes[alt_pos..alt_pos + len];
            let ref_x = ref_bytes[ref_pos];
            if let Some(p) = ins.iter().rposition(|&b| b == ref_x) {
                push_op(&mut result, p, 'I');
                push_op(&mut result, 1, '=');
                push_op(&mut result, len - p, 'I');
                ref_pos += 1;
                alt_pos += len + 1;
                i += 2;
                continue;
            }
        }

        // A single mismatch immediately followed by I(n): mirror of the above.
        if op == 'X' && len == 1 && i + 1 < ops.len() && ops[i + 1].1 == 'I' {
            let n = ops[i + 1].0;
            let ins = &alt_bytes[alt_pos + 1..alt_pos + 1 + n];
            let ref_x = ref_bytes[ref_pos];
            if let Some(p) = ins.iter().position(|&b| b == ref_x) {
                push_op(&mut result, 1 + p, 'I');
                push_op(&mut result, 1, '=');
                push_op(&mut result, n - p - 1, 'I');
                ref_pos += 1;
                alt_pos += 1 + n;
                i += 2;
                continue;
            }
        }

        // D(n) immediately followed by a single mismatch: mirror of I+X with
        // ref/alt roles swapped.
        if op == 'D' && i + 1 < ops.len() && ops[i + 1] == (1, 'X') {
            let del = &ref_bytes[ref_pos..ref_pos + len];
            let alt_x = alt_bytes[alt_pos];
            if let Some(q) = del.iter().rposition(|&b| b == alt_x) {
                push_op(&mut result, q, 'D');
                push_op(&mut result, 1, '=');
                push_op(&mut result, len - q, 'D');
                ref_pos += len + 1;
                alt_pos += 1;
                i += 2;
                continue;
            }
        }

        // A single mismatch immediately followed by D(n): mirror of X+I.
        if op == 'X' && len == 1 && i + 1 < ops.len() && ops[i + 1].1 == 'D' {
            let n = ops[i + 1].0;
            let del = &ref_bytes[ref_pos + 1..ref_pos + 1 + n];
            let alt_x = alt_bytes[alt_pos];
            if let Some(q) = del.iter().position(|&b| b == alt_x) {
                push_op(&mut result, 1 + q, 'D');
                push_op(&mut result, 1, '=');
                push_op(&mut result, n - q - 1, 'D');
                ref_pos += 1 + n;
                alt_pos += 1;
                i += 2;
                continue;
            }
        }

        match op {
            '=' | 'X' => {
                ref_pos += len;
                alt_pos += len;
            }
            'I' => alt_pos += len,
            'D' => ref_pos += len,
            _ => {}
        }
        push_op(&mut result, len, op);
        i += 1;
    }

    result
        .iter()
        .map(|(count, op)| format!("{}{}", count, op))
        .collect()
}

pub fn alignment_to_cigar(operations: &[AlignmentOperation]) -> String {
    let mut cigar: Vec<(usize, char)> = Vec::new();

    for op in operations {
        let cigar_op = match op {
            AlignmentOperation::Match => '=',
            AlignmentOperation::Subst => 'X',
            AlignmentOperation::Del => 'D',
            AlignmentOperation::Ins => 'I',
            AlignmentOperation::Xclip(_) => 'S',
            AlignmentOperation::Yclip(_) => 'S',
        };

        if !cigar.is_empty() && cigar.last().unwrap().1 == cigar_op {
            cigar.last_mut().unwrap().0 += 1;
        } else {
            cigar.push((1, cigar_op));
        }
    }

    cigar
        .iter()
        .map(|(count, op)| format!("{}{}", count, op))
        .collect()
}

pub fn gap_open_aligner(reference: &str, sequence: &str) -> String {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    // Create an aligner with the same scoring parameters
    let mut aligner = Aligner::with_capacity(sequence.len(), reference.len(), -5, -1, &score); // match_score=0, mismatch_score=-6, gap_open=-5, gap_extend=-3

    // Perform the alignment
    let alignment = aligner.global(sequence.as_bytes(), reference.as_bytes());
    let cigar = alignment_to_cigar(&alignment.operations);
    normalize_homopolymer_indels(&cigar, reference, sequence)
}

pub fn import_bed(bed_file: &String) -> Vec<(String, usize, usize)> {
    let mut bed_list = Vec::new();
    if bed_file.ends_with(".gz") {
        let file = File::open(bed_file).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line.unwrap();
            bed_list.push((
                line.split('\t').next().unwrap().to_string(),
                line.split('\t').nth(1).unwrap().parse().unwrap(),
                line.split('\t').nth(2).unwrap().parse().unwrap(),
            ));
        }
    } else {
        let file = File::open(bed_file).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line.unwrap();
            bed_list.push((
                line.split('\t').next().unwrap().to_string(),
                line.split('\t').nth(1).unwrap().parse().unwrap(),
                line.split('\t').nth(2).unwrap().parse().unwrap(),
            ));
        }
    }
    bed_list
}

pub fn process_cigar(cigar: &str) -> String {
    let mut out = String::new();
    let mut n = 0;

    for symbol in cigar.chars() {
        if symbol.is_ascii_digit() {
            n = 10 * n + symbol.to_digit(10).unwrap() as usize;
        } else {
            if n == 0 {
                out.push(symbol);
            } else {
                out.push_str(&symbol.to_string().repeat(n));
            }
            n = 0;
        }
    }

    out
}

pub fn combine_cigar(cigar: &str) -> String {
    if cigar.is_empty() {
        return String::new();
    }

    // Convert to Vec<char> for efficient indexing
    let mut chars: Vec<char> = cigar.chars().collect();
    chars.push('$'); // Add sentinel character to handle the last group

    let mut out = String::new();
    let mut start = 0;

    for i in 1..chars.len() {
        if chars[i - 1] != chars[i] {
            let length = i - start;
            out.push_str(&format!("{}{}", length, chars[i - 1]));
            start = i;
        }
    }

    out
}

pub fn jaccard_distance(vector1: &[bool], vector2: &[bool]) -> f64 {
    assert_eq!(
        vector1.len(),
        vector2.len(),
        "Vectors must have the same length"
    );

    let mut intersection_count = 0;
    let mut union_count = 0;

    for (a, b) in vector1.iter().zip(vector2.iter()) {
        if *a && *b {
            intersection_count += 1;
        }
        if *a || *b {
            union_count += 1;
        }
    }

    if union_count == 0 {
        return 0.0; // Both vectors are all zeros
    }

    1.0 - (intersection_count as f64 / union_count as f64)
}

/// Generate a null distribution through permutation testing
fn get_null_distribution(
    records: &Vec<String>,
    matrix: &Array2<f64>,
    permutation_round: usize,
) -> Vec<f64> {
    let bar = ProgressBar::new(permutation_round as u64);
    let summary_statistics = (0..permutation_round)
        .into_par_iter()
        .flat_map(|_| {
            bar.inc(1);
            let mut local_stats = Vec::new();
            let mut rng = SmallRng::from_os_rng();

            for (i, index) in records.iter().enumerate() {
                let vector = matrix.slice(s![i, ..]);

                // Create a shuffled copy of the vector
                let vector_data: Vec<f64> = vector.iter().copied().collect();
                let mut shuffled_data = vector_data.clone();
                shuffled_data.shuffle(&mut rng);
                let shuffled = Array1::from(shuffled_data);

                let mut all_coefficients = Vec::new();

                for (j, other_index) in records.iter().enumerate() {
                    if index == other_index {
                        continue;
                    }

                    let other_vector = matrix.slice(s![j, ..]);

                    let binary_vector: Vec<bool> = shuffled.iter().map(|&x| x > 0.5).collect();
                    let binary_other: Vec<bool> = other_vector.iter().map(|&x| x > 0.5).collect();

                    // Calculate Jaccard distance
                    let coor = (1.0 - jaccard_distance(&binary_vector, &binary_other)).abs();
                    all_coefficients.push(coor);
                }

                local_stats.push(all_coefficients.iter().sum());
            }
            local_stats
        })
        .collect::<Vec<f64>>();
    bar.finish();
    summary_statistics
}

/// Calculate statistics for observed data
fn calculate_observation_statistics(
    recordlist: &Vec<String>,
    index: usize,
    matrix: &Array2<f64>,
) -> f64 {
    let vector = &matrix.slice(s![index, ..]);
    let mut all_coefficients = Vec::new();

    for (i, other_index) in recordlist.iter().enumerate() {
        if i == index {
            continue;
        }

        let other_vector = &matrix.slice(s![i, ..]);

        // Convert arrays to binary vectors before calculating Jaccard distance
        let binary_vector: Vec<bool> = vector.iter().map(|&x| x > 0.5).collect();
        let binary_other: Vec<bool> = other_vector.iter().map(|&x| x > 0.5).collect();

        // Calculate Jaccard distance
        let coor = (1.0 - jaccard_distance(&binary_vector, &binary_other)).abs();
        all_coefficients.push(coor);
    }

    all_coefficients.iter().sum()
}

/// Calculate p-value using z-score approach
pub fn calculate_p_value(statistics: &[f64], observation: f64) -> f64 {
    let n = statistics.len() as f64;

    // Calculate mean
    let mu = statistics.iter().sum::<f64>() / n;

    // Calculate standard deviation
    let variance = statistics.iter().map(|&x| (x - mu).powi(2)).sum::<f64>() / n;
    let sigma = variance.sqrt();
    let z_score = (observation - mu) / sigma;

    // Calculate p-value using normal distribution CDF
    let normal = Normal::new(0.0, 1.0).unwrap();
    1.0 - normal.cdf(z_score)
}

pub fn permutation_test(
    matrix: &Array2<f64>,
    p_value_threshold: f64,
    permutation_round: usize,
    node_list: Vec<String>,
) -> Vec<String> {
    let bar = ProgressBar::new(node_list.len() as u64);
    let statistics = get_null_distribution(&node_list, matrix, permutation_round);
    let (indices, collected_values): (Vec<_>, Vec<_>) = (0..node_list.len())
        .into_par_iter()
        .map(|i| {
            bar.inc(1);
            let nodename = &node_list[i];
            let observation = calculate_observation_statistics(&node_list, i, matrix);
            let p_value = calculate_p_value(&statistics, observation);
            ((nodename.clone()), Some((p_value, nodename.clone())))
        })
        .unzip();

    let mut raw_p_values = Vec::new();
    let mut test_index = Vec::new();

    for item in collected_values.into_iter().flatten() {
        let (p_value, node_name) = item;
        if !p_value.is_nan() {
            raw_p_values.push(p_value);
            test_index.push(node_name);
        }
    }

    let mut excluded_index = Vec::new();

    if !raw_p_values.is_empty() {
        for (pi, p_value) in raw_p_values.iter().enumerate() {
            let test_index_value = &test_index[pi];
            if p_value > &p_value_threshold {
                excluded_index.push(test_index_value);
                debug!(
                    "excluded_index: {:?}, q_value: {:?}",
                    test_index_value, p_value
                );
            }
        }
    }

    bar.finish();

    let mut index_list = Vec::new();
    let mut f_node: Vec<String> = Vec::new();
    for (r, rindex) in node_list.iter().enumerate() {
        if !excluded_index.contains(&rindex) {
            index_list.push(r);
            f_node.push(rindex.clone());
        }
    }

    f_node
}

#[cfg(test)]
mod cigar_tests {
    use super::{gap_open_aligner, normalize_homopolymer_indels, parse_cigar};
    use std::collections::HashMap;

    #[test]
    fn parse_cigar_splits_runs() {
        assert_eq!(parse_cigar("3=2I1X"), vec![(3, '='), (2, 'I'), (1, 'X')]);
    }

    #[test]
    fn normalize_x_plus_i_homopolymer() {
        // Homopolymer T extension: aligner may emit 1X1I instead of 2I.
        let ref_seq = "ATTT";
        let alt_seq = "ATTTT";
        let raw = "3=1X1I";
        let normalized = normalize_homopolymer_indels(raw, ref_seq, alt_seq);
        assert_eq!(normalized, "3=1I1=");
        let (variants, _) = crate::call::get_variants_from_cigar(
            raw,
            "chr1",
            ref_seq,
            alt_seq,
            0,
            String::new(),
            &HashMap::new(),
            "n1",
        );
        assert!(
            variants.iter().all(|v| v.variant_type != "SNP"),
            "expected no spurious SNP, got {:?}",
            variants
        );
        let (variants_norm, _) = crate::call::get_variants_from_cigar(
            &normalized,
            "chr1",
            ref_seq,
            alt_seq,
            0,
            String::new(),
            &HashMap::new(),
            "n1",
        );
        assert_eq!(variants_norm.len(), 1);
        assert_eq!(variants_norm[0].variant_type, "INS");
    }

    #[test]
    fn normalize_i_plus_x_homopolymer() {
        let ref_seq = "ATTT";
        let alt_seq = "ATTTT";
        let raw = "3=1I1X";
        let normalized = normalize_homopolymer_indels(raw, ref_seq, alt_seq);
        assert_eq!(normalized, "4=1I");
    }

    #[test]
    fn normalize_leaves_real_snp_adjacent_to_indel() {
        // Genuine SNP next to indel should not be collapsed into the indel.
        let ref_seq = "AAAC";
        let alt_seq = "AAATAC";
        let raw = "3=1X1I1=";
        let normalized = normalize_homopolymer_indels(raw, ref_seq, alt_seq);
        assert_eq!(normalized, raw);
    }

    #[test]
    fn gap_open_aligner_avoids_snp_inside_indel() {
        let ref_seq = "ATTT";
        let alt_seq = "ATTTT";
        let cigar = gap_open_aligner(ref_seq, alt_seq);
        let (variants, _) = crate::call::get_variants_from_cigar(
            &cigar,
            "chr1",
            ref_seq,
            alt_seq,
            0,
            String::new(),
            &HashMap::new(),
            "n1",
        );
        assert!(
            variants.iter().all(|v| v.variant_type != "SNP"),
            "gap_open_aligner cigar {:?} produced SNP variants {:?}",
            cigar,
            variants
        );
    }
}
