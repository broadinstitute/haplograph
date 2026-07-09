use crate::intervals;
use crate::methyl;
use crate::util;
use anyhow::{Context, Result as AnyhowResult};
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use log::info;
use rust_htslib::bam;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{IndexedReader, Read as BamRead, Record as BamRecord};
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

/// Unique key for one BAM alignment record (qname + 0-based reference start).
pub fn alignment_bounds_key(qname: &str, align_pos: i64) -> String {
    format!("{qname}|{align_pos}")
}

pub fn normalized_read_slice_bounds(
    read_start: usize,
    read_end: usize,
    read_len: usize,
    _read_strand: &str,
) -> Option<(usize, usize)> {
    // BAM always stores sequences in forward-strand orientation for both plus and minus
    // strand reads (FLAG 0x10 means the stored SEQ is already the reverse-complement of
    // the original read, i.e. it is in the reference left→right direction).  Therefore
    // qpos always increases with reference position for both strands, start_pos < end_pos
    // in all normal cases, and no strand-specific coordinate transformation is needed.
    let pos_start = read_start.min(read_end);
    let pos_end = read_start.max(read_end);
    if pos_end == pos_start || pos_end > read_len {
        return None;
    }
    Some((pos_start, pos_end))
}

/// Per-alignment query bounds within a reference interval, plus window-edge alignment keys.
pub struct IntervalQueryBounds {
    /// First inclusive qpos in `[start, end)` keyed by `alignment_bounds_key`.
    pub per_align_start: HashMap<String, usize>,
    /// Last exclusive qpos in `[start, end)` keyed by `alignment_bounds_key`.
    pub per_align_end: HashMap<String, usize>,
    /// Alignment keys seen in the window with no valid qpos (fully deleted span).
    pub deletion_spanning: HashSet<String>,
    /// Read → alignment key active at the window start column.
    pub window_start_align: HashMap<String, String>,
    /// Read → alignment key active at the window end column.
    pub window_end_align: HashMap<String, String>,
}

/// One BAM fetch + pileup pass for an entire sub-locus; windows slice into this cache.
pub struct LocusPrefetch {
    pub chromosome: String,
    pub region_start: u64,
    pub region_end: u64,
    pub sampleid: String,
    records: Vec<BamRecord>,
    window_bounds: Vec<IntervalQueryBounds>,
}

struct BoundsBuilder {
    window_start: u64,
    window_end: u64,
    per_align_start: HashMap<String, usize>,
    per_align_end: HashMap<String, usize>,
    seen_in_window: HashSet<String>,
    window_start_align: HashMap<String, String>,
    window_end_align: HashMap<String, String>,
}

impl BoundsBuilder {
    fn new(window_start: u64, window_end: u64) -> Self {
        Self {
            window_start,
            window_end,
            per_align_start: HashMap::new(),
            per_align_end: HashMap::new(),
            seen_in_window: HashSet::new(),
            window_start_align: HashMap::new(),
            window_end_align: HashMap::new(),
        }
    }

    fn update_column(
        &mut self,
        ref_pos: u64,
        qname: &str,
        align_key: &str,
        qpos: Option<usize>,
        indel: bam::pileup::Indel,
    ) {
        self.seen_in_window.insert(align_key.to_string());
        if let Some(qpos) = qpos {
            self.per_align_start
                .entry(align_key.to_string())
                .or_insert(qpos);
            let excl = match indel {
                bam::pileup::Indel::Ins(len) => qpos + len as usize + 1,
                bam::pileup::Indel::Del(_) | bam::pileup::Indel::None => qpos + 1,
            };
            self.per_align_end.insert(align_key.to_string(), excl);
            if ref_pos == self.window_start {
                self.window_start_align
                    .entry(qname.to_string())
                    .or_insert_with(|| align_key.to_string());
            }
            self.window_end_align
                .insert(qname.to_string(), align_key.to_string());
        }
    }

    fn finish(self) -> IntervalQueryBounds {
        let deletion_spanning: HashSet<String> = self
            .seen_in_window
            .into_iter()
            .filter(|key| !self.per_align_start.contains_key(key))
            .collect();
        IntervalQueryBounds {
            per_align_start: self.per_align_start,
            per_align_end: self.per_align_end,
            deletion_spanning,
            window_start_align: self.window_start_align,
            window_end_align: self.window_end_align,
        }
    }
}

fn record_overlaps_interval(record: &BamRecord, start: u64, end: u64) -> bool {
    let ref_start = record.pos();
    let ref_end = record.reference_end();
    ref_start < end as i64 && ref_end > start as i64
}

/// Fetch all alignments and pileup bounds for every window in one indexed BAM pass.
pub fn build_locus_prefetch(
    bam: &mut IndexedReader,
    windows: &[(String, usize, usize)],
    sampleid: &String,
    primary_only: bool,
) -> AnyhowResult<LocusPrefetch> {
    anyhow::ensure!(
        !windows.is_empty(),
        "build_locus_prefetch requires at least one window"
    );

    let chromosome = windows[0].0.clone();
    let region_start = windows.first().unwrap().1 as u64;
    let region_end = windows.last().unwrap().2 as u64;

    bam.fetch((chromosome.as_bytes(), region_start, region_end))?;
    let mut records = Vec::new();
    for record_result in bam.records() {
        let record = record_result?;
        if primary_only && (record.is_secondary() || record.is_supplementary()) {
            continue;
        }
        records.push(record);
    }

    let mut builders: Vec<BoundsBuilder> = windows
        .iter()
        .map(|(_, ws, we)| BoundsBuilder::new(*ws as u64, *we as u64))
        .collect();

    bam.fetch((chromosome.as_bytes(), region_start, region_end))?;
    for pileup_result in bam.pileup() {
        let pileup = pileup_result?;
        let ref_pos = pileup.pos() as u64;
        if ref_pos < region_start || ref_pos >= region_end {
            continue;
        }
        let mut seen: HashSet<String> = HashSet::new();
        for alignment in pileup.alignments() {
            let record = alignment.record();
            let qname = String::from_utf8_lossy(record.qname()).into_owned();
            let align_key = alignment_bounds_key(&qname, record.pos());
            if seen.contains(&align_key) {
                continue;
            }
            seen.insert(align_key.clone());
            let qpos = alignment.qpos();
            let indel = alignment.indel();
            for (builder, (_, ws, we)) in builders.iter_mut().zip(windows.iter()) {
                if ref_pos >= *ws as u64 && ref_pos < *we as u64 {
                    builder.update_column(ref_pos, &qname, &align_key, qpos, indel);
                }
            }
        }
    }

    let window_bounds = builders.into_iter().map(|b| b.finish()).collect();

    info!(
        "Prefetched locus {}:{}-{} ({} records, {} windows)",
        chromosome,
        region_start,
        region_end,
        records.len(),
        windows.len()
    );

    Ok(LocusPrefetch {
        chromosome,
        region_start,
        region_end,
        sampleid: sampleid.clone(),
        records,
        window_bounds,
    })
}

/// Extract read sequences/coordinates for one window using a [`LocusPrefetch`] cache.
pub fn extract_window_coordinates_from_prefetch(
    prefetch: &LocusPrefetch,
    window_index: usize,
    window_start: u64,
    window_end: u64,
) -> AnyhowResult<(
    HashMap<String, (u64, u64)>,
    HashMap<String, Vec<u8>>,
    HashMap<String, Vec<u8>>,
    HashMap<String, String>,
    HashMap<String, HashMap<usize, f32>>,
)> {
    let bounds = prefetch
        .window_bounds
        .get(window_index)
        .with_context(|| format!("window index {window_index} out of range"))?;

    let chr = prefetch.chromosome.as_str();
    let sampleid = &prefetch.sampleid;

    let mut read_coordinates_formatted = HashMap::new();
    let mut read_sequence_dict_formatted = HashMap::new();
    let mut read_quality_dict_formatted = HashMap::new();
    let mut read_strand_dict_formatted = HashMap::new();
    let mut read_methyl_dict_formatted = HashMap::new();

    for record in &prefetch.records {
        if !record_overlaps_interval(record, window_start, window_end) {
            continue;
        }
        let read_name = String::from_utf8_lossy(record.qname()).into_owned();
        let align_key = alignment_bounds_key(&read_name, record.pos());

        if bounds.deletion_spanning.contains(&align_key) {
            let record_id = format!(
                "{}|{}",
                intervals::generate_read_name(&read_name, chr, window_start, window_end, sampleid),
                record.pos()
            );
            read_coordinates_formatted.insert(record_id.clone(), (0u64, 0u64));
            read_sequence_dict_formatted.insert(record_id.clone(), Vec::new());
            read_quality_dict_formatted.insert(record_id.clone(), Vec::new());
            read_strand_dict_formatted.insert(record_id.clone(), record.strand().to_string());
            read_methyl_dict_formatted.insert(record_id, HashMap::new());
            continue;
        }

        let Some(start_pos) = bounds.per_align_start.get(&align_key) else {
            continue;
        };
        let Some(end_pos) = bounds.per_align_end.get(&align_key) else {
            continue;
        };
        let read_seq = record.seq().as_bytes();
        let read_strand = record.strand().to_string();
        let read_len = read_seq.len();
        let Some((pos_start, pos_end)) =
            normalized_read_slice_bounds(*start_pos, *end_pos, read_len, &read_strand)
        else {
            continue;
        };

        let read_seq_final = String::from_utf8_lossy(&read_seq[pos_start..pos_end]).to_string();
        let record_id = format!(
            "{}|{}",
            intervals::generate_read_name(&read_name, chr, window_start, window_end, sampleid),
            record.pos()
        );
        read_coordinates_formatted.insert(record_id.clone(), (*start_pos as u64, *end_pos as u64));
        read_sequence_dict_formatted.insert(record_id.clone(), read_seq_final.clone().into_bytes());
        read_quality_dict_formatted.insert(record_id.clone(), record.qual().to_vec());
        read_strand_dict_formatted.insert(record_id.clone(), read_strand);
        read_methyl_dict_formatted.insert(
            record_id.clone(),
            methyl::get_methylation_read(record, *start_pos, *end_pos, 'm'),
        );
    }

    Ok((
        read_coordinates_formatted,
        read_quality_dict_formatted,
        read_sequence_dict_formatted,
        read_strand_dict_formatted,
        read_methyl_dict_formatted,
    ))
}

/// Build the in-window read sequence, concatenating across supplementary alignments when
/// the window start and end fall on different alignment records for the same read.
pub fn sequence_from_interval_bounds(
    records: &[BamRecord],
    qname: &str,
    bounds: &IntervalQueryBounds,
) -> Option<String> {
    let start_key = bounds.window_start_align.get(qname)?;
    let end_key = bounds.window_end_align.get(qname)?;

    if start_key == end_key {
        let record = records
            .iter()
            .find(|r| alignment_bounds_key(qname, r.pos()) == *start_key)?;
        let start_pos = bounds.per_align_start.get(start_key)?;
        let end_pos = bounds.per_align_end.get(end_key)?;
        let (pos_start, pos_end) =
            normalized_read_slice_bounds(*start_pos, *end_pos, record.seq().len(), "")?;
        return Some(
            String::from_utf8_lossy(&record.seq().as_bytes()[pos_start..pos_end]).to_string(),
        );
    }

    let start_ref_pos: i64 = start_key.rsplit('|').next()?.parse().ok()?;
    let end_ref_pos: i64 = end_key.rsplit('|').next()?.parse().ok()?;

    let mut segments: Vec<(i64, String)> = Vec::new();
    for record in records {
        if String::from_utf8_lossy(record.qname()) != qname {
            continue;
        }
        let align_key = alignment_bounds_key(qname, record.pos());
        let ref_pos = record.pos();
        if ref_pos < start_ref_pos || ref_pos > end_ref_pos {
            continue;
        }
        let Some(start_pos) = bounds.per_align_start.get(&align_key) else {
            continue;
        };
        let Some(end_pos) = bounds.per_align_end.get(&align_key) else {
            continue;
        };
        let Some((pos_start, pos_end)) =
            normalized_read_slice_bounds(*start_pos, *end_pos, record.seq().len(), "")
        else {
            continue;
        };
        let seq = String::from_utf8_lossy(&record.seq().as_bytes()[pos_start..pos_end]).to_string();
        segments.push((ref_pos, seq));
    }
    segments.sort_by_key(|(pos, _)| *pos);
    if segments.is_empty() {
        return None;
    }
    Some(segments.into_iter().map(|(_, seq)| seq).collect())
}

/// Single-pass pileup scan over `[start, end)`.
///
/// Returns `(start_dict, end_dict, deletion_spanning)` where:
///   - `start_dict`: each read's **first** valid inclusive qpos in the window.
///     Reads whose base at exactly `start` is deleted are still included using the
///     first non-deleted position they have within `[start, end)`.
///   - `end_dict`: each read's **last** valid exclusive qpos in the window.
///     Reads whose base at exactly `end-1` is deleted are still included using the
///     last non-deleted position they have within `[start, end)`.
///   - `deletion_spanning`: reads that appear in the pileup within `[start, end)` but
///     have `qpos = None` for every column — i.e. the entire query window falls inside
///     a large deletion of that read.  Callers should emit these reads with an empty
///     sequence (`""`).
///
/// Exclusive-end semantics per indel type at the last valid column:
///   - `None`   → `qpos + 1`         (include the matched base)
///   - `Ins(N)` → `qpos + N + 1`     (include anchor base + inserted bases)
///   - `Del(_)` → `qpos + 1`         (deletion starts outside the window; include the base)
pub fn record_interval_query_bounds_range(
    bam: &mut IndexedReader,
    start: u64,
    end: u64,
) -> AnyhowResult<IntervalQueryBounds> {
    let mut builder = BoundsBuilder::new(start, end);

    for p in bam.pileup() {
        let pileup = p?;
        let ref_pos = pileup.pos() as u64;
        if ref_pos < start {
            continue;
        }
        if ref_pos >= end {
            break;
        }
        let mut seen: HashSet<String> = HashSet::new();
        for alignment in pileup.alignments() {
            let record = alignment.record();
            let qname = String::from_utf8_lossy(record.qname()).into_owned();
            let align_key = alignment_bounds_key(&qname, record.pos());
            if seen.contains(&align_key) {
                continue;
            }
            seen.insert(align_key.clone());
            builder.update_column(
                ref_pos,
                &qname,
                &align_key,
                alignment.qpos(),
                alignment.indel(),
            );
        }
    }

    Ok(builder.finish())
}

// extract haplotypes from bam file
pub fn extract_haplotypes_from_bam(
    bam: &mut IndexedReader,
    chr: &str,
    start: u64,
    end: u64,
    sampleid: &String,
    primary_only: bool,
) -> AnyhowResult<HashMap<String, String>> {
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    // pb.set_message("Processing alignments...");

    let _ = bam.fetch((chr.as_bytes(), start, end));
    let mut records_by_qname: HashMap<String, Vec<BamRecord>> = HashMap::new();
    for record_result in bam.records() {
        let record = record_result?;
        if primary_only && (record.is_secondary() || record.is_supplementary()) {
            continue;
        }
        let read_name = String::from_utf8_lossy(record.qname()).into_owned();
        records_by_qname.entry(read_name).or_default().push(record);
    }

    let _ = bam.fetch((chr.as_bytes(), start, end));
    let bounds = record_interval_query_bounds_range(bam, start, end)?;

    let mut filtered_sequence_dict = HashMap::new();
    for (read_name, mut records) in records_by_qname {
        records.sort_by_key(|r| r.pos());
        let start_key = match bounds.window_start_align.get(&read_name) {
            Some(k) => k.clone(),
            None => continue,
        };
        if bounds.deletion_spanning.contains(&start_key) {
            let record_id = format!("{read_name}|{chr}:{start}-{end}|{sampleid}");
            filtered_sequence_dict.insert(record_id, String::new());
            continue;
        }

        let Some(read_seq_final) = sequence_from_interval_bounds(&records, &read_name, &bounds)
        else {
            continue;
        };

        let record_id = format!("{read_name}|{chr}:{start}-{end}|{sampleid}");
        filtered_sequence_dict.insert(record_id, read_seq_final);
        if filtered_sequence_dict.len() % 100 == 0 {
            pb.set_message(format!("Collected {} reads", filtered_sequence_dict.len()));
        }
    }

    pb.finish_with_message("Alignment processing completed");
    info!(
        "Extracted {} reads from region",
        filtered_sequence_dict.len()
    );

    Ok(filtered_sequence_dict)
}

// extract haplotypes from bam file
pub fn extract_haplotypes_coordinates_from_bam(
    bam: &mut IndexedReader,
    chr: &str,
    start: u64,
    end: u64,
    sampleid: &String,
    primary_only: bool,
) -> AnyhowResult<(
    HashMap<String, (u64, u64)>,
    HashMap<String, Vec<u8>>,
    HashMap<String, Vec<u8>>,
    HashMap<String, String>,
    HashMap<String, HashMap<usize, f32>>,
)> {
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    let _ = bam.fetch((chr.as_bytes(), start, end));
    let mut collected_records: Vec<BamRecord> = Vec::new();
    for record_result in bam.records() {
        let record = record_result?;
        if primary_only && (record.is_secondary() || record.is_supplementary()) {
            continue;
        }
        collected_records.push(record);
    }

    let _ = bam.fetch((chr.as_bytes(), start, end));
    let bounds = record_interval_query_bounds_range(bam, start, end)?;

    let mut read_coordinates_formatted = HashMap::new();
    let mut read_sequence_dict_formatted = HashMap::new();
    let mut read_quality_dict_formatted = HashMap::new();
    let mut read_strand_dict_formatted = HashMap::new();
    let mut read_methyl_dict_formatted = HashMap::new();
    for record in collected_records {
        let read_name = String::from_utf8_lossy(record.qname()).into_owned();
        let align_key = alignment_bounds_key(&read_name, record.pos());

        if bounds.deletion_spanning.contains(&align_key) {
            let record_id = format!(
                "{}|{}",
                intervals::generate_read_name(&read_name, chr, start, end, sampleid),
                record.pos()
            );
            read_coordinates_formatted.insert(record_id.clone(), (0u64, 0u64));
            read_sequence_dict_formatted.insert(record_id.clone(), Vec::new());
            read_quality_dict_formatted.insert(record_id.clone(), Vec::new());
            read_strand_dict_formatted.insert(record_id.clone(), record.strand().to_string());
            read_methyl_dict_formatted.insert(record_id, HashMap::new());
            continue;
        }

        let Some(start_pos) = bounds.per_align_start.get(&align_key) else {
            continue;
        };
        let Some(end_pos) = bounds.per_align_end.get(&align_key) else {
            continue;
        };
        let read_seq = record.seq().as_bytes();
        let read_strand = record.strand().to_string();
        let read_len = read_seq.len();
        let Some((pos_start, pos_end)) =
            normalized_read_slice_bounds(*start_pos, *end_pos, read_len, &read_strand)
        else {
            continue;
        };

        let read_seq_final = String::from_utf8_lossy(&read_seq[pos_start..pos_end]).to_string();

        let record_id = format!(
            "{}|{}",
            intervals::generate_read_name(&read_name, chr, start, end, sampleid),
            record.pos()
        );
        read_coordinates_formatted.insert(record_id.clone(), (*start_pos as u64, *end_pos as u64));
        read_sequence_dict_formatted.insert(record_id.clone(), read_seq_final.clone().into_bytes());
        read_quality_dict_formatted.insert(record_id.clone(), record.qual().to_vec());
        read_strand_dict_formatted.insert(record_id.clone(), read_strand);
        read_methyl_dict_formatted.insert(
            record_id.clone(),
            methyl::get_methylation_read(&record, *start_pos, *end_pos, 'm'),
        );

        if read_coordinates_formatted.len() % 100 == 0 {
            pb.set_message(format!(
                "Collected {} reads",
                read_coordinates_formatted.len()
            ));
        }
    }

    pb.finish_with_message("Alignment processing completed");
    info!(
        "Extracted {} reads from region",
        read_coordinates_formatted.len()
    );

    Ok((
        read_coordinates_formatted,
        read_quality_dict_formatted,
        read_sequence_dict_formatted,
        read_strand_dict_formatted,
        read_methyl_dict_formatted,
    ))
}

pub fn start(
    bam: &mut IndexedReader,
    chromosome: &str,
    start: usize,
    end: usize,
    primary_only: bool,
    output_path: String,
    sampleid: String,
    pileup: bool,
) -> AnyhowResult<()> {
    if pileup {
        let (
            reads,
            _read_coordinates,
            _read_sequence_dictionary,
            _read_methyl,
            _read_strand_dictionary,
        ) = intervals::extract_haplotypes_coordinates_from_bam_pileup(
            bam,
            chromosome,
            start as u64,
            end as u64,
            primary_only,
        )
        .unwrap();
        let outputfile = PathBuf::from(format!("{}.fasta", output_path));
        let mut file = File::create(outputfile)
            .with_context(|| format!("Failed to create output file: {}.fasta", output_path))?;

        for r in reads.iter() {
            let r_name = String::from_utf8_lossy(r.id().as_bytes()).to_string();
            let r_seq = String::from_utf8_lossy(r.seq()).to_string();
            writeln!(file, ">{}", r_name)?;
            // write the sequence in fasta format
            let seq_len = r_seq.len();
            let chars_per_line = 60;
            let full_lines = seq_len / chars_per_line;
            for i in 0..full_lines {
                let start = i * chars_per_line;
                let end = start + chars_per_line;
                writeln!(file, "{}", &r_seq[start..end])?;
            }
            if seq_len % chars_per_line != 0 {
                writeln!(file, "{}", &r_seq[full_lines * chars_per_line..])?;
            }
        }
    } else {
        let reads = extract_haplotypes_from_bam(
            bam,
            &chromosome,
            start as u64,
            end as u64,
            &sampleid,
            primary_only,
        )
        .unwrap();

        let output_p = PathBuf::from(format!("{}.fasta", output_path));
        if !reads.is_empty() {
            let _ = util::write_fasta(&reads, &output_p);
        }
    }

    Ok(())
}

#[cfg(test)]
mod extract_tests {
    use super::extract_haplotypes_from_bam;
    use rust_htslib::bam::IndexedReader;

    #[test]
    fn extract_haplotypes_from_bam_tp53_interval() {
        let bam_path = "HG00097_tp53_test.bam";
        if !std::path::Path::new(bam_path).exists() {
            return;
        }
        let mut bam = IndexedReader::from_path(bam_path).unwrap();
        let sampleid = "HG00097".to_string();
        let reads =
            extract_haplotypes_from_bam(&mut bam, "chr17", 7668752, 7668828, &sampleid, false)
                .unwrap();
        assert_eq!(reads.len(), 34);
    }

    #[test]
    fn prefetch_matches_per_window_extraction() {
        let bam_path = "HG00097_tp53_test.bam";
        if !std::path::Path::new(bam_path).exists() {
            return;
        }
        let mut bam = IndexedReader::from_path(bam_path).unwrap();
        let sampleid = "HG00097".to_string();
        let chr = "chr17";
        let start = 7668752_usize;
        let end = 7668828_usize;
        let windows = vec![(chr.to_string(), start, end)];

        let prefetch = super::build_locus_prefetch(&mut bam, &windows, &sampleid, false).unwrap();
        let (coords_p, _, seqs_p, _, _) =
            super::extract_window_coordinates_from_prefetch(&prefetch, 0, start as u64, end as u64)
                .unwrap();

        let mut bam2 = IndexedReader::from_path(bam_path).unwrap();
        let (coords_w, _, seqs_w, _, _) = super::extract_haplotypes_coordinates_from_bam(
            &mut bam2,
            chr,
            start as u64,
            end as u64,
            &sampleid,
            false,
        )
        .unwrap();

        assert_eq!(coords_p.len(), coords_w.len());
        assert_eq!(seqs_p.len(), seqs_w.len());
    }
}
