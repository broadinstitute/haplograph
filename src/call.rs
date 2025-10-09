
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
    pub node_id:String,
    pub haplotype_index: Option<Vec<usize>>,
}

pub fn get_variants_from_cigar(
    cigar: &str,
    ref_name: &str,
    ref_seq: &str,
    alt_seq: &str,
    ref_start: usize,
    allelecount: usize,
    node_id: &str,
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
                        node_id: node_id.to_string(),
                        haplotype_index: None,
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
                    node_id: node_id.to_string(),
                    haplotype_index: None,
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
                    node_id: node_id.to_string(),
                    haplotype_index: None,
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

    let mut collapsed_ad = HashMap::new();
    // let mut collapsed_af_node_id = HashMap::new();

    for current_var in variants {
        let chromosome = current_var.chromosome;
        let pos = current_var.pos;
        let ref_allele = current_var.ref_allele;
        let alt_allele = current_var.alt_allele;
        let variant_type = current_var.variant_type;
        let allele_count = current_var.allele_count;
        let node_id = current_var.node_id;
        let haplotype_index = current_var.haplotype_index;
        let key = (
            chromosome,
            pos,
            ref_allele.clone(),
            alt_allele.clone(),
            variant_type.clone(),
        );
        *collapsed_ad.entry(key).or_insert(0) += allele_count;
    }
    collapsed_ad
        .into_iter()
        .map(
            |((chromosome, pos, ref_allele, alt_allele, variant_type), allele_count)| Variant {
                chromosome,
                pos,
                ref_allele,
                alt_allele,
                variant_type,
                allele_count,
                node_id: "".to_string(),
                haplotype_index: None,
            },
        )
        .collect()
}


fn write_vcf(
    variants: &[Variant],
    coverage: &HashMap<usize, usize>,
    output_prefix: String,
    sample_id: &str,
    reference_seqs:&Vec<fastq::Record>,
    haplotype_number: usize,
    phase_variants: bool,
) -> AnyhowResult<()> {

    // Write VCF header
    // 1. Create a VCF header
    let mut header = bcf::Header::new();
    for record in reference_seqs.iter() {
        let referencename = record.id().to_string();
        let referencelength = record.seq().len() as u64;
        header.push_record(format!("##contig=<ID={},length={}>\n", referencename, referencelength).as_bytes());
    }

    header.push_record(format!("##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">\n").as_bytes());
    header.push_record(format!("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n").as_bytes());
    header.push_record(format!("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n").as_bytes());
    header.push_record(format!("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Alternative Allele Depth\">\n").as_bytes());
    header.push_record(format!("##FORMAT=<ID=VAF,Number=1,Type=Float,Description=\"Variant Allele Frequency\">\n").as_bytes());
    header.push_sample(sample_id.as_bytes());

    // Write VCF
    // 2. Open a compressed VCF writer
    // Create a BCF writer for the compressed VCF.
    let output_file = PathBuf::from(format!("{}.vcf.gz", output_prefix));
    let mut writer = bcf::Writer::from_path(output_file.clone(), &header, false, bcf::Format::Vcf).expect("Failed to create BCF writer");

    // Sort variants by chromosome and position
    let mut sorted_variants = variants.to_vec();
    sorted_variants.sort_by(|a, b| {
        // First compare chromosomes
        a.chromosome
            .cmp(&b.chromosome)
            // Then compare positions
            .then(a.pos.cmp(&b.pos))
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
    //sort records_by_pos by the key
    let mut record_by_pos_keys = records_by_pos.keys().collect::<Vec<_>>();
    record_by_pos_keys.sort_by(|a, b| {
        a.1.cmp(&b.1)
        .then(a.2.cmp(&b.2))
    });
    // println!("records_by_pos: {:?}", record_by_pos_keys);

    
    // vcf_records
    for key in record_by_pos_keys.iter() {
        let mut var_list = records_by_pos.get(key).unwrap().clone();
        var_list.sort_by(|a, b| a.alt_allele.cmp(&b.alt_allele));
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
        record.set_pos(variant_pos as i64 -1);
        record.set_id(b".");
        let mut all_variant_alleles = vec![variant_ref_allele] ;
        all_variant_alleles.extend(variant_alt_allele);
        record.set_alleles(&all_variant_alleles.iter().map(|v| v.as_bytes()).collect::<Vec<&[u8]>>()).expect("Failed to set alleles");
        record.push_format_integer(b"DP", &[*coverage.get(&variant_pos).unwrap_or(&0) as i32]).expect("Failed to set DP format field");
        record.push_format_integer(b"AD", &var_list.iter().map(|v| v.allele_count as i32).collect::<Vec<i32>>()).expect("Failed to set AD format field");
        record.push_format_float(b"VAF", &allele_frequency).expect("Failed to set VAF format field");
        let mut genotype_list = Vec::new();
        // phase variants
        if phase_variants {
            let mut phase = false;
            for haplotype_index in 1..=haplotype_number {
                let mut found = false;
                for (index, variant) in var_list.iter().enumerate() {
                    let haplotype_index_list = variant.clone().haplotype_index.unwrap();
                    // println!("haplotype_index_: {}", haplotype_index_);
                    for hap_index in haplotype_index_list.iter() {
                        if *hap_index == haplotype_index {
                            if haplotype_index == 1 {
                                genotype_list.push(GenotypeAllele::Unphased(index as i32 + 1));
                            } else {
                                genotype_list.push(GenotypeAllele::Phased(index as i32 + 1));
                            }
                            found = true;
                            phase = true;
                            break;
                        } 
                    }
                }
                if !found {
                    if phase {
                        genotype_list.push(GenotypeAllele::Phased(0));
                    } else {
                        genotype_list.push(GenotypeAllele::Unphased(0));
                    }
                }
            }
            if !phase {
                record.push_info_flag(b"SOMATIC").expect("Failed to set SOMATIC info field");
            }
        } else {
            if var_list.len() == 1 {
                genotype_list.push(GenotypeAllele::Unphased(0));
            } 
            for (index, variant) in var_list.iter().enumerate() {
                genotype_list.push(GenotypeAllele::Unphased(index as i32 + 1));
            }
            
        }

        record.push_genotypes(&genotype_list.as_slice()).expect("Failed to set genotype field");


        writer.write(&record).expect("Failed to write record");
    }

    
    info!("Created compressed VCF: {}", &output_file.display());
    
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

/// Phase variants and write phased VCF
pub fn Phase_germline_variants(
    graph_filename: &PathBuf, 
    Variants: &Vec<Variant>,
    haplotype_number: usize,
) -> AnyhowResult<Vec<Variant>> {
    let (node_info, edge_info) = asm::load_graph(graph_filename).unwrap();
    // establish the phasing information
    // let all_sequences = asm::traverse_graph(&node_info, &edge_info, true, haplotype_number).unwrap();
    let (germline_nodes_sorted, germline_edge_info, node_haplotype, haplotype_reads) = asm::find_germline_nodes(&node_info, &edge_info, haplotype_number);

    let mut collapsed_variants = HashMap::new();
    for variant in Variants.iter() {
        let key = (variant.chromosome.clone(), variant.pos, variant.ref_allele.clone(), variant.alt_allele.clone(), variant.variant_type.clone());
        collapsed_variants.entry(key).or_insert(Vec::new()).push((variant.allele_count.clone(), variant.node_id.clone()));
    }

    let mut phased_variants = Vec::new();
    for (key_t,node_list) in collapsed_variants.iter() {
        // println!("variant: {:?}", variant.node_id);
        let chromosome = key_t.0.clone();
        let pos = key_t.1.clone();
        let ref_allele = key_t.2.clone();
        let alt_allele = key_t.3.clone();
        let variant_type = key_t.4.clone();
        let allele_count = node_list.iter().map(|x| x.0).sum();
        let node_id_list = node_list.iter().map(|x| x.1.clone()).collect::<Vec<String>>();
        let mut haplotype_index = Vec::new();
        for node_id in node_id_list.iter() {
            if node_haplotype.contains_key(node_id) {
                // println!("node_id: {:?}, haplotype_map: {:?}", node_id, node_haplotype.get(node_id).unwrap());
                haplotype_index.extend(node_haplotype.get(node_id).unwrap().iter().map(|x| x + 1));
            }
        }
        // println!("haplotype_index: {:?}", haplotype_index);

        phased_variants.push(Variant {
            chromosome: chromosome,
            pos: pos,
            ref_allele: ref_allele,
            alt_allele: alt_allele,
            variant_type: variant_type,
            allele_count: allele_count,
            node_id: node_id_list.join(","),
            haplotype_index: Some(haplotype_index),
        });
    }
    Ok(phased_variants)
}



pub fn start(graph_filename: &PathBuf, reference_seqs: &Vec<fastq::Record>, sampleid: &String, output_prefix: &String, haplotype_number: usize, phase_variants: bool) -> AnyhowResult<()> {
    let (node_info, _) = asm::load_graph(graph_filename).unwrap();
    let mut all_variants = Vec::new();
    let mut coverage = HashMap::new();
    for (node_id, n_info) in node_info.iter() {
        let cigar = n_info.cigar.clone();
        let alt_seq = n_info.seq.clone();
        let support_reads = n_info.support_reads.clone();
        // transform support_reads to usize, drop the quote before and after the number
        // let support_reads = support_reads.trim_matches('"').parse::<usize>().unwrap();
        // let support_reads = support_reads.parse::<usize>().unwrap();
        // println!("support_reads: {}", support_reads);
        let parts: Vec<&str> = node_id.split(".").collect();
        let locus = parts[1];
        let (chromosome, start, end) = util::split_locus(locus.to_string());
        let full_ref_seq_ = reference_seqs.iter().find(|r| r.id().to_string() == chromosome).unwrap().seq();
        let full_ref_seq = String::from_utf8_lossy(full_ref_seq_).to_string();
        let ref_seq = full_ref_seq[start..end].to_string();
        let (variants, poscounts) = get_variants_from_cigar(cigar.as_str(), &chromosome, ref_seq.as_str(), alt_seq.as_str(), start, support_reads, &node_id);
        all_variants.extend(variants);

        for (pos, count) in poscounts.iter() {
            *coverage.entry(*pos).or_insert(0) += count;
        }

    }

    let variants = if phase_variants {
        Phase_germline_variants(&graph_filename, &all_variants, haplotype_number)?
    } else {
        collapse_identical_records(all_variants)
    };
    
    // Write VCF file    
    write_vcf(
        &variants,
        &coverage,
        output_prefix.clone(),
        sampleid,
        reference_seqs,
        haplotype_number,
        phase_variants
    )?; 
    info!("Total variants: {}", variants.len());

    Ok(())
}