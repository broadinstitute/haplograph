use crate::asm;
use crate::util;
use anyhow::{bail, Context, Result as AnyhowResult};
use bio::io::fasta::Reader as FastaReader;
use log::{info, warn};
use rust_htslib::bcf::{self, Read};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::ops::Deref;
use std::path::{Path, PathBuf};

/// Default overlap between adjacent maximal-locus windows (matches `main.rs`).
pub const DEFAULT_OVERLAP_BP: usize = 500;

#[derive(Debug, Clone)]
struct ParsedHaplotype {
    hap_index: usize,
    chromosome: String,
    ref_start: usize,
    ref_end: usize,
    sequence: String,
    path: Vec<String>,
    read_names: HashSet<String>,
}

#[derive(Debug, Clone)]
struct Segment {
    index: usize,
    prefix: String,
    ref_start: usize,
    ref_end: usize,
    haplotypes: HashMap<usize, ParsedHaplotype>,
}

/// One genomic chunk in merge order; `segment` is `None` when assembly was skipped or failed.
#[derive(Debug, Clone)]
struct MergeChunk {
    chrom: String,
    ref_start: usize,
    ref_end: usize,
    segment: Option<Segment>,
}

#[derive(Debug, Clone)]
struct MethylRow {
    chrom: String,
    ref_start: usize,
    ref_end: usize,
    mod_rate: f32,
    asm_start: usize,
    asm_end: usize,
    motif: String,
    coverage: usize,
}

/// Plan overlapping sub-loci the same way as `haplograph haplograph` in `main.rs`.
pub fn plan_locus_chunks(
    chromosome: &str,
    locus_start: usize,
    locus_end: usize,
    maximal_locus_size: usize,
    overlap_bp: usize,
) -> AnyhowResult<Vec<(String, usize, usize)>> {
    let locus_len = locus_end.saturating_sub(locus_start);
    if locus_len == 0 {
        bail!("Empty locus interval");
    }
    if locus_len <= maximal_locus_size {
        return Ok(vec![(chromosome.to_string(), locus_start, locus_end)]);
    }
    if maximal_locus_size <= overlap_bp {
        bail!(
            "maximal_locus_size ({maximal_locus_size}) must be greater than overlap ({overlap_bp} bp)"
        );
    }
    let step = maximal_locus_size - overlap_bp;
    let mut chunks = Vec::new();
    let mut start_pos = locus_start;
    while start_pos < locus_end {
        let end_pos = (start_pos + maximal_locus_size).min(locus_end);
        chunks.push((chromosome.to_string(), start_pos, end_pos));
        if end_pos >= locus_end {
            break;
        }
        start_pos += step;
    }
    Ok(chunks)
}

/// Discover `{output_prefix}_0`, `{output_prefix}_1`, … segment prefixes in order.
pub fn discover_segment_prefixes(output_prefix: &str) -> AnyhowResult<Vec<String>> {
    let mut prefixes = Vec::new();
    let mut index = 0usize;
    loop {
        let candidate = format!("{output_prefix}_{index}");
        let fasta_path = PathBuf::from(format!("{candidate}.fasta"));
        if fasta_path.exists() {
            prefixes.push(candidate);
            index += 1;
        } else if index == 0 {
            bail!(
                "No segment FASTA found at {}.fasta — run haplograph on a large locus first",
                candidate
            );
        } else {
            break;
        }
    }
    Ok(prefixes)
}

fn parse_fasta_header(header: &str) -> AnyhowResult<(String, usize, usize, usize, Vec<String>)> {
    let parts: Vec<&str> = header.split('\t').collect();
    let id_part = parts.first().context("Missing FASTA header id")?;
    let path = if parts.len() >= 3 {
        parts[2].split('|').map(|s| s.to_string()).collect()
    } else {
        Vec::new()
    };

    let mut id_fields = id_part.split('.');
    let coord_part = id_fields
        .next()
        .context("Missing coordinate field in FASTA header")?;
    let hap_index = id_fields
        .next()
        .context("Missing haplotype index in FASTA header")?
        .parse::<usize>()?;

    let mut coord_split = coord_part.split(':');
    let chromosome = coord_split
        .next()
        .context("Missing chromosome in FASTA header")?
        .to_string();
    let span = coord_split.next().context("Missing span in FASTA header")?;
    let mut span_split = span.split('-');
    let ref_start = span_split
        .next()
        .context("Missing start in FASTA header")?
        .parse::<usize>()?;
    let ref_end = span_split
        .next()
        .context("Missing end in FASTA header")?
        .parse::<usize>()?;

    Ok((chromosome, ref_start, ref_end, hap_index, path))
}

fn path_read_names(node_info: &HashMap<String, asm::NodeInfo>, path: &[String]) -> HashSet<String> {
    let mut reads = HashSet::new();
    for node in path {
        if node_info.contains_key(node) {
            reads.extend(asm::get_read_name_list(node_info, node.clone()));
        }
    }
    reads
}

fn node_overlaps_interval(node_id: &str, overlap_start: usize, overlap_end: usize) -> bool {
    let locus = match node_id.split('.').nth(1) {
        Some(value) => value,
        None => return false,
    };
    let (_, node_start, node_end) = util::split_locus(locus.to_string());
    node_start < overlap_end && node_end > overlap_start
}

fn hap_overlap_reads(
    hap: &ParsedHaplotype,
    node_info: &HashMap<String, asm::NodeInfo>,
    overlap_start: usize,
    overlap_end: usize,
) -> HashSet<String> {
    let mut reads = HashSet::new();
    for node in &hap.path {
        if node_overlaps_interval(node, overlap_start, overlap_end) {
            reads.extend(asm::get_read_name_list(node_info, node.clone()));
        }
    }
    reads
}

fn load_chromosome_reference(reference_fa: &str, chromosome: &str) -> AnyhowResult<String> {
    let (_, chr_seqs) = util::get_ref_seq_from_chromosome(&reference_fa.to_string(), chromosome);
    let ref_record = chr_seqs
        .first()
        .with_context(|| format!("Chromosome {chromosome} not found in reference"))?;
    Ok(String::from_utf8_lossy(ref_record.seq()).to_string())
}

fn reference_slice(full_ref: &str, start: usize, end: usize) -> String {
    let end = end.min(full_ref.len());
    let start = start.min(end);
    full_ref[start..end].to_string()
}

fn try_load_segment(prefix: &str, index: usize) -> Option<Segment> {
    let fasta_path = PathBuf::from(format!("{prefix}.fasta"));
    if !fasta_path.exists() {
        return None;
    }
    match load_segment(prefix, index) {
        Ok(segment) => Some(segment),
        Err(err) => {
            warn!("Failed to load segment {prefix}: {err}");
            None
        }
    }
}

fn discover_all_segments(output_prefix: &str, max_scan: usize) -> Vec<Segment> {
    let mut segments = Vec::new();
    for index in 0..max_scan {
        let prefix = format!("{output_prefix}_{index}");
        if let Some(segment) = try_load_segment(&prefix, index) {
            segments.push(segment);
        }
    }
    segments
}

fn build_merge_plan(output_prefix: &str, planned: &[(String, usize, usize)]) -> Vec<MergeChunk> {
    let max_scan = planned.len().saturating_add(32);
    let loaded_by_coords: HashMap<(usize, usize), Segment> =
        discover_all_segments(output_prefix, max_scan)
            .into_iter()
            .map(|segment| ((segment.ref_start, segment.ref_end), segment))
            .collect();

    planned
        .iter()
        .map(|(chrom, start, end)| {
            let segment = loaded_by_coords.get(&(*start, *end)).cloned();
            if segment.is_none() {
                warn!(
                    "No assembly for {chrom}:{start}-{end}; both haplotypes will use reference in merged FASTA"
                );
            }
            MergeChunk {
                chrom: chrom.clone(),
                ref_start: *start,
                ref_end: *end,
                segment,
            }
        })
        .collect()
}

fn load_segment(prefix: &str, index: usize) -> AnyhowResult<Segment> {
    let fasta_path = PathBuf::from(format!("{prefix}.fasta"));
    let gfa_path = PathBuf::from(format!("{prefix}.gfa"));
    let node_info = if gfa_path.exists() {
        asm::load_graph(&gfa_path)
            .with_context(|| format!("Failed to load GFA for segment {prefix}"))?
            .0
    } else {
        warn!("Missing GFA for segment {prefix}; read-based phasing will be limited");
        HashMap::new()
    };

    let reader = FastaReader::from_file(&fasta_path)
        .with_context(|| format!("Failed to open {fasta_path:?}"))?;
    let mut haplotypes = HashMap::new();
    let mut ref_start = usize::MAX;
    let mut ref_end = 0usize;

    for record in reader.records() {
        let record = record.with_context(|| format!("Failed to read FASTA from {fasta_path:?}"))?;
        let header = record.id();
        let (chromosome, start, end, hap_index, path) = parse_fasta_header(header)?;
        ref_start = ref_start.min(start);
        ref_end = ref_end.max(end);
        let read_names = path_read_names(&node_info, &path);
        haplotypes.insert(
            hap_index,
            ParsedHaplotype {
                hap_index,
                chromosome,
                ref_start: start,
                ref_end: end,
                sequence: String::from_utf8_lossy(record.seq()).to_string(),
                path,
                read_names,
            },
        );
    }

    if haplotypes.is_empty() {
        bail!("No haplotype records found in {fasta_path:?}");
    }

    Ok(Segment {
        index,
        prefix: prefix.to_string(),
        ref_start,
        ref_end,
        haplotypes,
    })
}

/// Map global haplotype indices (from segment 0) to local hap indices in `right`.
/// Uses read overlap in the shared reference interval between adjacent segments.
fn match_haplotypes_at_overlap(
    left: &Segment,
    right: &Segment,
    left_gfa: &HashMap<String, asm::NodeInfo>,
    right_gfa: &HashMap<String, asm::NodeInfo>,
    number_of_haplotypes: usize,
) -> AnyhowResult<Vec<usize>> {
    let overlap_start = right.ref_start;
    let overlap_end = left.ref_end.min(right.ref_end);
    if overlap_start >= overlap_end {
        warn!(
            "Segments {} and {} do not share reference overlap; using identity phasing",
            left.prefix, right.prefix
        );
        return Ok((0..number_of_haplotypes).collect());
    }

    let left_reads: HashMap<usize, HashSet<String>> = left
        .haplotypes
        .iter()
        .map(|(hap, parsed)| {
            (
                *hap,
                hap_overlap_reads(parsed, left_gfa, overlap_start, overlap_end),
            )
        })
        .collect();
    let right_reads: HashMap<usize, HashSet<String>> = right
        .haplotypes
        .iter()
        .map(|(hap, parsed)| {
            (
                *hap,
                hap_overlap_reads(parsed, right_gfa, overlap_start, overlap_end),
            )
        })
        .collect();

    let empty = HashSet::new();
    let mut best_mapping = (0..number_of_haplotypes).collect::<Vec<_>>();
    let mut best_score = 0usize;

    let permutations: Vec<Vec<usize>> = if number_of_haplotypes == 1 {
        vec![vec![0]]
    } else if number_of_haplotypes == 2 {
        vec![vec![0, 1], vec![1, 0]]
    } else {
        bail!("merge currently supports 1 or 2 haplotypes, got {number_of_haplotypes}");
    };

    for perm in permutations {
        let score: usize = (0..number_of_haplotypes)
            .map(|global_hap| {
                left_reads
                    .get(&global_hap)
                    .unwrap_or(&empty)
                    .intersection(right_reads.get(&perm[global_hap]).unwrap_or(&empty))
                    .count()
            })
            .sum();
        if score > best_score {
            best_score = score;
            best_mapping = perm;
        }
    }

    info!(
        "Phasing {} -> {} at overlap {}:{}-{} (read-overlap score={})",
        left.prefix,
        right.prefix,
        left.haplotypes
            .values()
            .next()
            .map(|h| h.chromosome.as_str())
            .unwrap_or("?"),
        overlap_start,
        overlap_end,
        best_score
    );
    Ok(best_mapping)
}

fn ref_coord_to_asm_offset(hap: &ParsedHaplotype, ref_pos: usize) -> usize {
    if ref_pos <= hap.ref_start {
        return 0;
    }
    if ref_pos >= hap.ref_end {
        return hap.sequence.len();
    }
    let ref_span = hap.ref_end.saturating_sub(hap.ref_start).max(1);
    ((ref_pos - hap.ref_start) * hap.sequence.len()) / ref_span
}

fn stitch_haplotype_pair(
    left: &ParsedHaplotype,
    right: &ParsedHaplotype,
    overlap_start: usize,
    overlap_end: usize,
) -> String {
    let trim_left = ref_coord_to_asm_offset(left, overlap_start);
    let skip_right = ref_coord_to_asm_offset(right, overlap_end);
    let mut stitched = left.sequence[..trim_left.min(left.sequence.len())].to_string();
    if skip_right < right.sequence.len() {
        stitched.push_str(&right.sequence[skip_right..]);
    }
    stitched
}

fn build_global_to_local_maps(
    merge_chunks: &[MergeChunk],
    number_of_haplotypes: usize,
) -> AnyhowResult<Vec<Vec<usize>>> {
    let mut global_to_local = vec![(0..number_of_haplotypes).collect::<Vec<_>>()];
    for idx in 1..merge_chunks.len() {
        let left = merge_chunks[idx - 1].segment.as_ref();
        let right = merge_chunks[idx].segment.as_ref();
        let mapping = match (left, right) {
            (Some(left_seg), Some(right_seg))
                if !left_seg.prefix.is_empty() && !right_seg.prefix.is_empty() =>
            {
                let left_gfa_path = PathBuf::from(format!("{}.gfa", left_seg.prefix));
                let right_gfa_path = PathBuf::from(format!("{}.gfa", right_seg.prefix));
                if left_gfa_path.exists() && right_gfa_path.exists() {
                    let left_gfa = asm::load_graph(&left_gfa_path)?.0;
                    let right_gfa = asm::load_graph(&right_gfa_path)?.0;
                    match_haplotypes_at_overlap(
                        left_seg,
                        right_seg,
                        &left_gfa,
                        &right_gfa,
                        number_of_haplotypes,
                    )?
                } else {
                    (0..number_of_haplotypes).collect()
                }
            }
            _ => (0..number_of_haplotypes).collect(),
        };
        global_to_local.push(mapping);
    }
    Ok(global_to_local)
}

fn haplotype_for_chunk(
    chunk: &MergeChunk,
    global_hap: usize,
    local_hap: usize,
    full_ref: &str,
) -> ParsedHaplotype {
    if let Some(segment) = &chunk.segment {
        if let Some(hap) = segment.haplotypes.get(&local_hap) {
            return hap.clone();
        }
        warn!(
            "Missing haplotype {local_hap} in {}; using reference for haplotype {global_hap}",
            segment.prefix
        );
    }
    ParsedHaplotype {
        hap_index: local_hap,
        chromosome: chunk.chrom.clone(),
        ref_start: chunk.ref_start,
        ref_end: chunk.ref_end,
        sequence: reference_slice(full_ref, chunk.ref_start, chunk.ref_end),
        path: Vec::new(),
        read_names: HashSet::new(),
    }
}

fn append_reference_gap(stitched: &mut String, full_ref: &str, gap_start: usize, gap_end: usize) {
    if gap_start < gap_end {
        stitched.push_str(&reference_slice(full_ref, gap_start, gap_end));
    }
}

fn stitch_fasta(
    merge_chunks: &[MergeChunk],
    global_to_local: &[Vec<usize>],
    number_of_haplotypes: usize,
    full_ref: &str,
    locus_start: usize,
    locus_end: usize,
) -> AnyhowResult<HashMap<String, String>> {
    if merge_chunks.is_empty() {
        bail!("No merge chunks to stitch");
    }

    let mut output = HashMap::new();
    for global_hap in 0..number_of_haplotypes {
        let mut stitched = String::new();
        let mut tracked_ref_start = locus_start;
        let mut tracked_ref_end = locus_start;

        if merge_chunks[0].ref_start > locus_start {
            append_reference_gap(
                &mut stitched,
                full_ref,
                locus_start,
                merge_chunks[0].ref_start,
            );
        }

        for (chunk_idx, chunk) in merge_chunks.iter().enumerate() {
            let local_hap = global_to_local[chunk_idx][global_hap];
            let hap = haplotype_for_chunk(chunk, global_hap, local_hap, full_ref);

            if stitched.is_empty() {
                stitched = hap.sequence.clone();
                tracked_ref_start = chunk.ref_start;
                tracked_ref_end = chunk.ref_end;
                continue;
            }

            if chunk_idx > 0 {
                let prev_end = merge_chunks[chunk_idx - 1].ref_end;
                if chunk.ref_start > prev_end {
                    append_reference_gap(&mut stitched, full_ref, prev_end, chunk.ref_start);
                }
            }

            let overlap_start = chunk.ref_start;
            let overlap_end = merge_chunks[chunk_idx - 1].ref_end.min(chunk.ref_end);
            let left_hap = ParsedHaplotype {
                hap_index: global_hap,
                chromosome: chunk.chrom.clone(),
                ref_start: tracked_ref_start,
                ref_end: tracked_ref_end,
                sequence: stitched.clone(),
                path: Vec::new(),
                read_names: HashSet::new(),
            };
            stitched = stitch_haplotype_pair(&left_hap, &hap, overlap_start, overlap_end);
            tracked_ref_end = chunk.ref_end;
        }

        if let Some(last_chunk) = merge_chunks.last() {
            if last_chunk.ref_end < locus_end {
                append_reference_gap(&mut stitched, full_ref, last_chunk.ref_end, locus_end);
            }
        }

        let header = format!(
            "{}:{}-{}.{global_hap}",
            merge_chunks[0].chrom, locus_start, locus_end
        );
        output.insert(header, stitched);
    }
    Ok(output)
}

fn load_methyl_bed(path: &Path) -> AnyhowResult<Vec<MethylRow>> {
    if !path.exists() {
        return Ok(Vec::new());
    }
    let file = File::open(path).with_context(|| format!("Failed to open methyl bed {path:?}"))?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            continue;
        }
        rows.push(MethylRow {
            chrom: fields[0].to_string(),
            ref_start: fields[1].parse()?,
            ref_end: fields[2].parse()?,
            mod_rate: fields[3].parse()?,
            asm_start: fields[4].parse()?,
            asm_end: fields[5].parse()?,
            motif: fields[6].to_string(),
            coverage: fields[7].parse()?,
        });
    }
    Ok(rows)
}

fn write_methyl_bed(path: &Path, haplotype_index: usize, rows: &[MethylRow]) -> AnyhowResult<()> {
    let mut file = File::create(path)?;
    writeln!(file, "##fileformat=BED")?;
    writeln!(file, "##haplotype={}", haplotype_index + 1)?;
    writeln!(
        file,
        "#CHROM\tRef_start\tRef_end\tMod_rate\tAsm_start\tAsm_end\tMotif\tCoverage"
    )?;
    for row in rows {
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.chrom,
            row.ref_start,
            row.ref_end,
            row.mod_rate,
            row.asm_start,
            row.asm_end,
            row.motif,
            row.coverage
        )?;
    }
    Ok(())
}

fn merge_methyl_beds(
    merge_chunks: &[MergeChunk],
    global_to_local: &[Vec<usize>],
    stitched_sequences: &HashMap<String, String>,
    number_of_haplotypes: usize,
    output_prefix: &str,
) -> AnyhowResult<()> {
    for global_hap in 0..number_of_haplotypes {
        let mut merged_rows = Vec::new();
        let mut asm_offset = 0usize;
        let mut prev_assembly_end: Option<usize> = None;
        for (chunk_idx, chunk) in merge_chunks.iter().enumerate() {
            let segment = match &chunk.segment {
                Some(segment) => segment,
                None => continue,
            };
            let local_hap = global_to_local[chunk_idx][global_hap];
            let bed_path = PathBuf::from(format!("{}.Hap.{}.bed", segment.prefix, local_hap));
            let mut rows = load_methyl_bed(&bed_path)?;
            let ref_min = prev_assembly_end.unwrap_or(chunk.ref_start);
            rows.retain(|row| row.ref_start >= ref_min);
            for row in rows.iter_mut() {
                row.asm_start += asm_offset;
                row.asm_end += asm_offset;
            }
            if chunk_idx + 1 < merge_chunks.len() {
                let overlap_end = chunk.ref_end;
                let overlap_start = merge_chunks[chunk_idx + 1].ref_start;
                let unique_len = overlap_end.saturating_sub(overlap_start);
                asm_offset += unique_len;
            } else if let Some(seq) = stitched_sequences.values().next() {
                asm_offset = seq.len();
            }
            merged_rows.extend(rows);
            prev_assembly_end = Some(chunk.ref_end);
        }
        merged_rows.sort_by_key(|row| row.ref_start);
        let bed_path = PathBuf::from(format!("{output_prefix}.Hap.{global_hap}.bed"));
        write_methyl_bed(&bed_path, global_hap, &merged_rows)?;
        info!(
            "Wrote merged methylation BED for haplotype {} ({} sites)",
            global_hap,
            merged_rows.len()
        );
    }
    Ok(())
}

fn swap_genotype_haplotypes(record: &mut bcf::Record, number_of_haplotypes: usize) {
    if number_of_haplotypes < 2 {
        return;
    }
    let genotypes = match record.genotypes() {
        Ok(g) => g,
        Err(_) => return,
    };
    let gt = genotypes.get(0);
    let mut alleles = gt.deref().clone();
    if alleles.len() >= 2 {
        alleles.swap(0, 1);
        let _ = record.push_genotypes(&alleles);
    }
}

fn merge_vcf_files(
    merge_chunks: &[MergeChunk],
    global_to_local: &[Vec<usize>],
    number_of_haplotypes: usize,
    sample_id: &str,
    output_prefix: &str,
    vcf_kind: &str,
) -> AnyhowResult<()> {
    let template_prefix = merge_chunks
        .iter()
        .find_map(|chunk| chunk.segment.as_ref().map(|segment| segment.prefix.clone()))
        .with_context(|| format!("No {vcf_kind} template segment found for merge"))?;
    let template_path = format!("{template_prefix}.{vcf_kind}.vcf.gz");
    let template_reader = bcf::Reader::from_path(&template_path)
        .with_context(|| format!("Failed to open template VCF {template_path}"))?;
    let mut header = bcf::Header::from_template(template_reader.header());
    header.push_sample(sample_id.as_bytes());

    let output_path = format!("{output_prefix}.{vcf_kind}.vcf.gz");
    let mut writer = bcf::Writer::from_path(&output_path, &header, false, bcf::Format::Vcf)
        .with_context(|| format!("Failed to create merged VCF {output_path}"))?;

    let mut seen_positions: HashSet<(String, i64, String)> = HashSet::new();
    let mut prev_assembly_end: Option<usize> = None;

    for (chunk_idx, chunk) in merge_chunks.iter().enumerate() {
        let segment = match &chunk.segment {
            Some(segment) => segment,
            None => continue,
        };
        let vcf_path = format!("{}.{}.vcf.gz", segment.prefix, vcf_kind);
        if !Path::new(&vcf_path).exists() {
            warn!("Missing {vcf_path}; skipping in merge");
            continue;
        }
        let mut reader = bcf::Reader::from_path(&vcf_path)
            .with_context(|| format!("Failed to open {vcf_path}"))?;
        let local_to_global: HashMap<usize, usize> = global_to_local[chunk_idx]
            .iter()
            .enumerate()
            .map(|(global_hap, local_hap)| (*local_hap, global_hap))
            .collect();
        let flip = chunk_idx > 0
            && number_of_haplotypes == 2
            && local_to_global.get(&0) == Some(&1)
            && local_to_global.get(&1) == Some(&0);

        let ref_min = prev_assembly_end.unwrap_or(chunk.ref_start);

        for record_result in reader.records() {
            let mut record = record_result?;
            let chrom = match record.rid() {
                Some(rid) => String::from_utf8_lossy(record.header().rid2name(rid).unwrap_or(b"?"))
                    .to_string(),
                None => "?".to_string(),
            };
            let pos = record.pos() + 1;
            if (pos as usize) < ref_min {
                continue;
            }
            let ref_allele = record
                .alleles()
                .first()
                .map(|allele| String::from_utf8_lossy(allele).to_string())
                .unwrap_or_default();
            let key = (chrom.clone(), pos, ref_allele.clone());
            if seen_positions.contains(&key) {
                continue;
            }
            seen_positions.insert(key);

            if flip {
                swap_genotype_haplotypes(&mut record, number_of_haplotypes);
            }
            writer.write(&record)?;
        }
        prev_assembly_end = Some(chunk.ref_end);
    }

    info!("Wrote merged {vcf_kind} VCF: {output_path}");
    Ok(())
}

/// Stitch per-window haplograph outputs into full-length FASTA, VCF, and methylation BED files.
///
/// Expects segment outputs named `{output_prefix}_0`, `{output_prefix}_1`, … each with
/// `.fasta`, `.gfa`, `.germline.vcf.gz`, `.somatic.vcf.gz`, and `.Hap.N.bed` files.
pub fn start(
    output_prefix: &str,
    locus: &str,
    number_of_haplotypes: usize,
    sample_id: &str,
    reference_fa: &str,
    maximal_locus_size: usize,
    overlap_bp: usize,
) -> AnyhowResult<()> {
    let (chromosome, locus_start, locus_end) = util::split_locus(locus.to_string());
    let planned = plan_locus_chunks(
        &chromosome,
        locus_start,
        locus_end,
        maximal_locus_size,
        overlap_bp,
    )?;
    info!(
        "Merging {} planned sub-loci for locus {}:{}-{} (overlap={} bp)",
        planned.len(),
        chromosome,
        locus_start,
        locus_end,
        overlap_bp
    );

    let merge_chunks = build_merge_plan(output_prefix, &planned);
    let assembled_count = merge_chunks
        .iter()
        .filter(|chunk| chunk.segment.is_some())
        .count();
    info!(
        "Found assembly for {assembled_count}/{} planned sub-loci",
        merge_chunks.len()
    );
    if assembled_count == 0 {
        bail!("No segment assemblies found under prefix {output_prefix}_*");
    }

    let full_ref = load_chromosome_reference(reference_fa, &chromosome)?;
    let global_to_local = build_global_to_local_maps(&merge_chunks, number_of_haplotypes)?;
    let stitched_haplotypes = stitch_fasta(
        &merge_chunks,
        &global_to_local,
        number_of_haplotypes,
        &full_ref,
        locus_start,
        locus_end,
    )?;
    let fasta_path = PathBuf::from(format!("{output_prefix}.fasta"));
    util::write_fasta(&stitched_haplotypes, &fasta_path)?;
    info!(
        "Wrote merged haplotype FASTA: {} ({} records)",
        fasta_path.display(),
        stitched_haplotypes.len()
    );

    merge_methyl_beds(
        &merge_chunks,
        &global_to_local,
        &stitched_haplotypes,
        number_of_haplotypes,
        output_prefix,
    )?;

    merge_vcf_files(
        &merge_chunks,
        &global_to_local,
        number_of_haplotypes,
        sample_id,
        output_prefix,
        "germline",
    )?;
    merge_vcf_files(
        &merge_chunks,
        &global_to_local,
        number_of_haplotypes,
        sample_id,
        output_prefix,
        "somatic",
    )?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_fasta_header_extracts_coords_and_path() {
        let header = "chr6:1000-2000.1\tSupports:42\tH.chr6:1000-1100.0|H.chr6:1100-1200.0";
        let (chrom, start, end, hap, path) = parse_fasta_header(header).unwrap();
        assert_eq!(chrom, "chr6");
        assert_eq!(start, 1000);
        assert_eq!(end, 2000);
        assert_eq!(hap, 1);
        assert_eq!(path.len(), 2);
    }

    #[test]
    fn stitch_haplotype_pair_trims_overlap() {
        let left = ParsedHaplotype {
            hap_index: 0,
            chromosome: "chr1".to_string(),
            ref_start: 1000,
            ref_end: 2000,
            sequence: "A".repeat(1000),
            path: Vec::new(),
            read_names: HashSet::new(),
        };
        let right = ParsedHaplotype {
            hap_index: 0,
            chromosome: "chr1".to_string(),
            ref_start: 1500,
            ref_end: 2500,
            sequence: "C".repeat(1000),
            path: Vec::new(),
            read_names: HashSet::new(),
        };
        let stitched = stitch_haplotype_pair(&left, &right, 1500, 2000);
        assert_eq!(stitched.len(), 1000);
        assert_eq!(stitched.matches('A').count(), 500);
        assert_eq!(stitched.matches('C').count(), 500);
    }

    #[test]
    fn stitch_fasta_uses_reference_for_missing_chunk() {
        let full_ref = "A".repeat(1000);
        let merge_chunks = vec![
            MergeChunk {
                chrom: "chr1".to_string(),
                ref_start: 100,
                ref_end: 300,
                segment: Some(Segment {
                    index: 0,
                    prefix: "s0".to_string(),
                    ref_start: 100,
                    ref_end: 300,
                    haplotypes: HashMap::from([(
                        0,
                        ParsedHaplotype {
                            hap_index: 0,
                            chromosome: "chr1".to_string(),
                            ref_start: 100,
                            ref_end: 300,
                            sequence: "G".repeat(200),
                            path: Vec::new(),
                            read_names: HashSet::new(),
                        },
                    )]),
                }),
            },
            MergeChunk {
                chrom: "chr1".to_string(),
                ref_start: 500,
                ref_end: 700,
                segment: None,
            },
        ];
        let global_to_local = vec![vec![0], vec![0]];
        let stitched =
            stitch_fasta(&merge_chunks, &global_to_local, 1, &full_ref, 100, 700).unwrap();
        let seq = stitched.values().next().unwrap();
        // 200 bp assembly + 200 bp gap ref + 200 bp missing-chunk ref
        assert_eq!(seq.len(), 600);
        assert!(seq.starts_with('G'));
        assert!(seq.ends_with('A'));
    }

    #[test]
    fn match_haplotypes_prefers_read_overlap() {
        let left = Segment {
            index: 0,
            prefix: "left".to_string(),
            ref_start: 1000,
            ref_end: 2000,
            haplotypes: HashMap::from([
                (
                    0,
                    ParsedHaplotype {
                        hap_index: 0,
                        chromosome: "chr1".to_string(),
                        ref_start: 1000,
                        ref_end: 2000,
                        sequence: String::new(),
                        path: vec!["H.chr1:1500-1600.0".to_string()],
                        read_names: HashSet::from(["r1".to_string(), "r2".to_string()]),
                    },
                ),
                (
                    1,
                    ParsedHaplotype {
                        hap_index: 1,
                        chromosome: "chr1".to_string(),
                        ref_start: 1000,
                        ref_end: 2000,
                        sequence: String::new(),
                        path: vec!["H.chr1:1500-1600.1".to_string()],
                        read_names: HashSet::from(["r3".to_string()]),
                    },
                ),
            ]),
        };
        let right = Segment {
            index: 1,
            prefix: "right".to_string(),
            ref_start: 1500,
            ref_end: 2500,
            haplotypes: HashMap::from([
                (
                    0,
                    ParsedHaplotype {
                        hap_index: 0,
                        chromosome: "chr1".to_string(),
                        ref_start: 1500,
                        ref_end: 2500,
                        sequence: String::new(),
                        path: vec!["H.chr1:1500-1600.0".to_string()],
                        read_names: HashSet::from(["r9".to_string()]),
                    },
                ),
                (
                    1,
                    ParsedHaplotype {
                        hap_index: 1,
                        chromosome: "chr1".to_string(),
                        ref_start: 1500,
                        ref_end: 2500,
                        sequence: String::new(),
                        path: vec!["H.chr1:1500-1600.1".to_string()],
                        read_names: HashSet::from(["r1".to_string(), "r2".to_string()]),
                    },
                ),
            ]),
        };

        let left_gfa = HashMap::from([
            (
                "H.chr1:1500-1600.0".to_string(),
                asm::NodeInfo {
                    seq: String::new(),
                    cigar: String::new(),
                    support_reads: 2,
                    allele_frequency: "0.5".to_string(),
                    read_names: "r1,r2".to_string(),
                    methyl_info: HashMap::new(),
                },
            ),
            (
                "H.chr1:1500-1600.1".to_string(),
                asm::NodeInfo {
                    seq: String::new(),
                    cigar: String::new(),
                    support_reads: 1,
                    allele_frequency: "0.5".to_string(),
                    read_names: "r3".to_string(),
                    methyl_info: HashMap::new(),
                },
            ),
        ]);
        let right_gfa = HashMap::from([
            (
                "H.chr1:1500-1600.0".to_string(),
                asm::NodeInfo {
                    seq: String::new(),
                    cigar: String::new(),
                    support_reads: 1,
                    allele_frequency: "0.5".to_string(),
                    read_names: "r9".to_string(),
                    methyl_info: HashMap::new(),
                },
            ),
            (
                "H.chr1:1500-1600.1".to_string(),
                asm::NodeInfo {
                    seq: String::new(),
                    cigar: String::new(),
                    support_reads: 2,
                    allele_frequency: "0.5".to_string(),
                    read_names: "r1,r2".to_string(),
                    methyl_info: HashMap::new(),
                },
            ),
        ]);

        let mapping = match_haplotypes_at_overlap(&left, &right, &left_gfa, &right_gfa, 2).unwrap();
        assert_eq!(mapping, vec![1, 0]);
    }
}
