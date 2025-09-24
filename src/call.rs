
use log::info;
use anyhow::{Result as AnyhowResult, Context};
use bio::io::fastq;
use std::path::PathBuf;
use std::collections::HashMap;
use crate::util;
use crate::asm;
use rust_htslib::bcf::{self, index as BcfIndex, record::GenotypeAllele};




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
                        pos: pos + 1,
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
) -> AnyhowResult<()> {

    // Write VCF header
    // 1. Create a VCF header
    let mut header = bcf::Header::new();
    for record in reference_seqs.iter() {
        let referencename = record.id().to_string();
        let referencelength = record.seq().len() as u64;
        header.push_record(format!("##contig=<ID={},length={}>\n", referencename, referencelength).as_bytes());
    }

    header.push_record(format!("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n").as_bytes());
    header.push_record(format!("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n").as_bytes());
    header.push_record(format!("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n").as_bytes());
    header.push_record(format!("##FORMAT=<ID=VAF,Number=1,Type=Float,Description=\"Variant Allele Frequency\">\n").as_bytes());
    header.push_sample(sample_id.as_bytes());

    // Write VCF
    // 2. Open a compressed VCF writer
    // Create a BCF writer for the compressed VCF.
    let mut writer = bcf::Writer::from_path(output_file, &header, false, bcf::Format::Vcf).expect("Failed to create BCF writer");

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

    // Merge multi-allelic variants
    let mut records_by_pos = HashMap::new();
    for var in sorted_variants {
        let chrom = var.chromosome.clone();
        let pos = var.pos.clone();
        let ref_allele = var.ref_allele.clone();
        let key = (chrom, pos, ref_allele);

        records_by_pos.entry(key)
            .or_insert_with(Vec::new)
            .push(var.clone());
    }
    
    // vcf_records
    for (key, var_list) in records_by_pos.iter() {
        let variant_chromosome = key.0.clone();
        let variant_pos = key.1.clone();
        let variant_ref_allele = key.2.clone();
        let variant_alt_allele = var_list.iter().map(|v| v.alt_allele.clone()).collect::<Vec<String>>();
        let read_depth = coverage.get(&variant_pos).unwrap_or(&0);
        let allele_frequency = if *read_depth == 0 {
            vec![0.0; var_list.len()]
        } else {
            var_list.iter().map(|v| v.allele_count as f32 / *read_depth as f32).collect::<Vec<f32>>()
        };
    
        // Create a variant record.
        let mut record = writer.empty_record();
        let header_view = writer.header();
        let rid = header_view.name2rid(variant_chromosome.as_bytes()).unwrap();        
        record.set_rid(Some(rid));
        record.set_pos(variant_pos as i64);
        record.set_id(b".");
        let mut all_variant_alleles = vec![variant_ref_allele] ;
        all_variant_alleles.extend(variant_alt_allele);
        record.set_alleles(&all_variant_alleles.iter().map(|v| v.as_bytes()).collect::<Vec<&[u8]>>()).expect("Failed to set alleles");
        record.push_format_integer(b"DP", &[*coverage.get(&variant_pos).unwrap_or(&0) as i32]).expect("Failed to set DP format field");
        record.push_format_integer(b"AD", &var_list.iter().map(|v| v.allele_count as i32).collect::<Vec<i32>>()).expect("Failed to set AD format field");
        record.push_format_float(b"VAF", &allele_frequency).expect("Failed to set VAF format field");
        let mut genotype_list = Vec::new();
        if allele_frequency.len() == 1 {
            if allele_frequency[0] > 0.99 {
                genotype_list.push(GenotypeAllele::Unphased(1));
                genotype_list.push(GenotypeAllele::Unphased(1));
            } else {
                genotype_list.push(GenotypeAllele::Unphased(0));
                genotype_list.push(GenotypeAllele::Unphased(1));
            }
        } else {
            for af in 0..allele_frequency.len()+1 {
                genotype_list.push(GenotypeAllele::Unphased(af as i32));
                
            }

        }
        record.push_genotypes(&genotype_list.as_slice()).expect("Failed to set genotype field");


        writer.write(&record).expect("Failed to write record");
    }

    
    info!("Created compressed VCF: {}", output_file.display());
    
    // Build tabix index for the compressed VCF file
    let idx = std::ptr::null();
    println!("Building index file: ");
    let rs = unsafe {
        rust_htslib::htslib::bcf_index_build3(
            rust_htslib::utils::path_to_cstring(&output_file).unwrap().as_ptr(),
            idx,
            0,
            4 as i32,
    )};
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