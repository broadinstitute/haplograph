use log::{info, warn};
use anyhow::{Result as AnyhowResult, Context};
use bio::io::fastq;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, record::Aux};
use std::path::{PathBuf};
use std::collections::{HashMap, BTreeMap, HashSet};
use std::fs::File;
use crate::util;
use crate::asm;
use rust_htslib::bcf::{self, Reader as BcfReader, Writer as BcfWriter, index as BcfIndex};
use flate2::write::GzEncoder;
use flate2::Compression;

use std::io::{Read, Write};
use std::process::Command;



#[derive(Debug, Clone)]
pub struct Variant {
    pub chromosome: String,
    pub pos: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub variant_type: String,
    pub allele_count: usize,
}

pub fn get_variants_from_cigar(
    cigar: &str,
    ref_name: &str,
    ref_seq: &str,
    alt_seq: &str,
    ref_start: usize,
    allelecount: usize,
) -> (Vec<Variant>, HashMap<usize, usize>) {
    let mut poscount = HashMap::new();
    let mut variants = Vec::new();
    let mut ref_pos = 0;
    let mut alt_pos = 0;

    let mut operations: Vec<(usize, char)> = Vec::new();
    let mut num = String::new();

    for c in cigar.trim_matches('"').chars() {
        // trim_matches removes quotes at start and end
        if c.is_digit(10) {
            num.push(c);
        } else {
            let number = num.parse::<usize>().expect("number {}");
            operations.push((number, c));
            num.clear();
        }
    }

    for (length, op) in operations {
        match op {
            '=' => {
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    *poscount.entry(pos).or_insert(0) += allelecount;
                }
                ref_pos += length;
                alt_pos += length;
            }
            'X' => {
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    *poscount.entry(pos).or_insert(0) += allelecount;
                    let ref_allele = &ref_seq[ref_pos + i..ref_pos + i + 1];
                    let alt_allele = &alt_seq[alt_pos + i..alt_pos + i + 1];
                    variants.push(Variant {
                        chromosome: ref_name.to_string(),
                        pos: pos,
                        ref_allele: ref_allele.to_string(),
                        alt_allele: alt_allele.to_string(),
                        variant_type: "SNP".to_string(),
                        allele_count: allelecount,
                    });
                }
                ref_pos += length;
                alt_pos += length;
            }
            'I' => {
                let pos = ref_start + ref_pos;
                // *poscount.entry(pos).or_insert(0) += allelecount;
                let ref_allele = if ref_pos > 0 {
                    match ref_seq.get(ref_pos - 1..ref_pos) {
                        Some(allele) => allele,
                        None => {
                            println!("{} {} {}", ref_seq, ref_seq.len(), ref_pos);
                            "-"
                        }
                    }
                } else {
                    "-"
                };

                let alt_allele = match alt_seq.get(alt_pos - 1..alt_pos + length) {
                    Some(allele) => allele,
                    None => {
                        println!("{} {} {}", alt_seq, alt_seq.len(), alt_pos);
                        "-"
                    }
                };
                variants.push(Variant {
                    chromosome: ref_name.to_string(),
                    pos: pos,
                    ref_allele: ref_allele.to_string(),
                    alt_allele: alt_allele.to_string(),
                    variant_type: "INS".to_string(),
                    allele_count: allelecount,
                });
                alt_pos += length;
            }

            'D' => {
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    *poscount.entry(pos).or_insert(0) += allelecount;
                }
                let pos = ref_start + ref_pos;
                
                let ref_allele = if ref_pos > 0 { match ref_seq.get(ref_pos - 1..ref_pos + length) {
                    Some(allele) => allele,
                    None => {
                        println!("{} {} {}", ref_seq, ref_seq.len(), ref_pos);
                        "-"
                    }
                } } else {
                   &format!("N{}", match ref_seq.get(ref_pos..ref_pos + length) {
                        Some(allele) => allele,
                        None => {
                            println!("{} {} {}", ref_seq, ref_seq.len(), ref_pos);
                            "-"
                        }
                    })};

                let alt_allele = if alt_pos > 0 {
                    match alt_seq.get(alt_pos - 1..alt_pos) {
                        Some(allele) => allele,
                        None => {
                            println!("{} {} {}", alt_seq, alt_seq.len(), alt_pos);
                            "-"
                        }
                    }
                } else {
                    "-"
                };

                if ref_pos > 0 && alt_pos > 0 {
                    if let (Some(r), Some(a)) = (
                        ref_seq.get(ref_pos - 1..ref_pos),
                        alt_seq.get(alt_pos - 1..alt_pos),
                    ) {
                        if r != a {
                            println!("{} {} {}", r, a, cigar);
                        }
                    }
                }

                variants.push(Variant {
                    chromosome: ref_name.to_string(),
                    pos: pos,
                    ref_allele: ref_allele.to_string(),
                    alt_allele: alt_allele.to_string(),
                    variant_type: "DEL".to_string(),
                    allele_count: allelecount,
                });
                ref_pos += length;
            }
            _ => (), //others skip
        }
    }
    (variants, poscount)
}

fn collapse_identical_records(variants: Vec<Variant>) -> Vec<Variant> {
    if variants.is_empty() {
        return Vec::new();
    }

    let mut collapsed = HashMap::new();

    for current_var in variants {
        let chromosome = current_var.chromosome;
        let pos = current_var.pos;
        let ref_allele = current_var.ref_allele;
        let alt_allele = current_var.alt_allele;
        let variant_type = current_var.variant_type;
        let allele_count = current_var.allele_count;

        let key = (
            chromosome,
            pos,
            ref_allele.clone(),
            alt_allele.clone(),
            variant_type.clone(),
        );
        *collapsed.entry(key).or_insert(0) += allele_count;
    }
    collapsed
        .into_iter()
        .map(
            |((chromosome, pos, ref_allele, alt_allele, variant_type), allele_count)| Variant {
                chromosome,
                pos,
                ref_allele,
                alt_allele,
                variant_type,
                allele_count,
            },
        )
        .collect()
}

fn format_vcf_record(variant: &Variant, coverage: HashMap<usize, usize>) -> String {
    // Add AC (allele count) to INFO field
    let read_depth = coverage.get(&variant.pos).unwrap_or(&0);
    let allele_frequency = if *read_depth == 0 {
        0.0
    } else {
        variant.allele_count as f32 / *read_depth as f32
    };

    let info = format!("DP={}", read_depth);
    let format: String = format!("GT:AD:VAF");
    // should be modified based on the phasing information
    let genotype: String = format!("1/1");
    let sample: String = format!("{}:{}:{}", genotype, variant.allele_count, allele_frequency);
    
    
    match variant.variant_type.as_str() {
        "SNP" => format!(
            "{}\t{}\t.\t{}\t{}\t.\t.\t{}\t{}\t{}",
            variant.chromosome,
            variant.pos + 1,
            variant.ref_allele,
            variant.alt_allele,
            info,
            format,
            sample
        ),
        "INS" => format!(
            "{}\t{}\t.\t{}\t{}\t.\t.\t{}\t{}\t{}",
            variant.chromosome,
            variant.pos, variant.ref_allele, variant.alt_allele, info, format, sample
        ),
        "DEL" => format!(
            "{}\t{}\t.\t{}\t{}\t.\t.\t{}\t{}\t{}",
            variant.chromosome,
            variant.pos, variant.ref_allele, variant.alt_allele, info, format, sample
        ),
        _ => panic!("Unknown variant type"),
    }
}

fn write_vcf(
    variants: &[Variant],
    coverage: &HashMap<usize, usize>,
    output_file: &PathBuf,
    sample_id: &str,
    reference_seqs:&Vec<fastq::Record>,
) -> std::io::Result<()> {
    // let mut file = File::create(Path::new(output_file))?;
   // Create compressed VCF file
   let compressed_filename = output_file.with_extension("vcf.gz");
   let compressed_file = File::create(&compressed_filename)?;
   let mut gz_encoder = GzEncoder::new(compressed_file, Compression::default());
//    gz_encoder.write_all(vcf_content.as_bytes())?;
//    gz_encoder.finish()?;

    // Write VCF header
    gz_encoder.write_all(b"##fileformat=VCFv4.2\n")?;
    
    for record in reference_seqs.iter() {
        let referencename = record.id().to_string();
        let referencelength = record.seq().len();
        gz_encoder.write_all(format!("##contig=<ID={},length={}>\n", referencename, referencelength).as_bytes())?;
    }
    
    gz_encoder.write_all(b"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")?;
    gz_encoder.write_all(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")?;
    gz_encoder.write_all(b"##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n")?;
    gz_encoder.write_all(b"##FORMAT=<ID=VAF,Number=1,Type=Float,Description=\"Variant Allele Frequency\">\n")?;
    gz_encoder.write_all(format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n", sample_id).as_bytes())?;

    // Sort variants by chromosome and position
    let mut sorted_variants = variants.to_vec();
    sorted_variants.sort_by(|a, b| {
        // First compare chromosomes
        a.chromosome
            .cmp(&b.chromosome)
            // Then compare positions
            .then(a.pos.cmp(&b.pos))
            // Then compare variant types (to ensure consistent ordering)
            .then(a.variant_type.cmp(&b.variant_type))
            // Then compare ref alleles
            .then(a.ref_allele.cmp(&b.ref_allele))
            // Then compare alt alleles
            .then(a.alt_allele.cmp(&b.alt_allele))
    });

    // Write variant records
    for variant in sorted_variants {
        gz_encoder.write_all(format!("{}\n", format_vcf_record(&variant, coverage.clone())).as_bytes())?;
    }
    gz_encoder.finish()?;

    
    info!("Created compressed VCF: {}", compressed_filename.display());
    // info!("Created index file: {}", index_filename.display());

    Ok(())
}


pub fn start(graph_filename: &PathBuf, reference_seqs: &Vec<fastq::Record>, sampleid: &String, output_prefix: &PathBuf) -> AnyhowResult<()> {
    let (node_info, _) = asm::load_graph(graph_filename).unwrap();
    let mut all_variants = Vec::new();
    let mut coverage = HashMap::new();
    for (node_id, n_info) in node_info.iter() {
        let cigar = n_info.cigar.clone();
        let alt_seq = n_info.seq.clone();
        let support_reads = n_info.support_reads.clone();
        // transform support_reads to usize, drop the quote before and after the number
        let support_reads = support_reads.trim_matches('"').parse::<usize>().unwrap();
        // let support_reads = support_reads.parse::<usize>().unwrap();
        // println!("support_reads: {}", support_reads);
        let parts: Vec<&str> = node_id.split(".").collect();
        let locus = parts[1];
        let (chromosome, start, end) = util::split_locus(locus.to_string());
        let full_ref_seq_ = reference_seqs.iter().find(|r| r.id().to_string() == chromosome).unwrap().seq();
        let full_ref_seq = String::from_utf8_lossy(full_ref_seq_).to_string();
        let ref_seq = full_ref_seq[start..end].to_string();
        let (variants, poscounts) = get_variants_from_cigar(cigar.as_str(), &chromosome, ref_seq.as_str(), alt_seq.as_str(), start, support_reads);
        all_variants.extend(variants);

        for (pos, count) in poscounts.iter() {
            *coverage.entry(*pos).or_insert(0) += count;
        }

    }
    let collapsed_variants = collapse_identical_records(all_variants);
    
    // Write VCF file    
    write_vcf(
        &collapsed_variants,
        &coverage,
        &output_prefix,
        sampleid,
        reference_seqs
    )?;
    info!("Total variants: {}", collapsed_variants.len());

    Ok(())
}