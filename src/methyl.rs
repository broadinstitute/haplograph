use rust_htslib::bam::{self, IndexedReader, Read, Reader, Record, Read as BamRead, record::Aux, pileup::Pileup};
use std::collections::{HashMap, HashSet};
use log::{warn};

pub fn get_methylation_read(r: &Record, start: usize, end: usize, mod_char: char) -> HashMap<usize, f32> {
    let mut methyl_pos_dict: HashMap<usize, f32> = HashMap::new();
    let forward_sequence = String::from_utf8_lossy(&r.seq().as_bytes()).to_string();
    
    // Check for modification data
    if let Ok(mods) = r.basemods_iter() {
        // Iterate over the modification types
        for res in mods {
            if let Ok( (position, m) ) = res {
                if m.modified_base as u8 as char != mod_char{
                    continue
                }
                // let strand = mod_metadata.strand;
                if position < start as i32 || position > end as i32 {
                    continue
                }
                // println!("position: {}, start: {}, end: {}", position, start, end);
                let pos_usize = position as usize;
                let qual = m.qual as f32 / 255.0;
                let motif: String = if r.is_reverse(){
                    forward_sequence[pos_usize-1..pos_usize +1 ].to_string()
                    
                }else{
                    forward_sequence[pos_usize..pos_usize +2 ].to_string()
                };
                    
                if motif == "CG" {
                    if r.is_reverse(){
                        if pos_usize as i32 - start as i32 - 1 > 0 {
                            methyl_pos_dict.insert(pos_usize - 1 - start as usize, qual);
                        }
                        
                    }else{
                        methyl_pos_dict.insert(pos_usize - start as usize, qual);

                    }
                    
                }else{
                    warn!("Unexpected motif: {},{},{}", pos_usize, qual, motif);
                    continue
                }

            }
        }                    

    }
    
    methyl_pos_dict
}

pub fn aggregate_methylation_reads(methyl_all_reads: HashMap<String,HashMap<usize, f32>>, min_methyl_threshold: f32) -> HashMap<usize, f32> {
    let mut aggregated_methyl_pos_dict: HashMap<usize, Vec<f32>> = HashMap::new();
    for (read_name, methyl_info) in methyl_all_reads.iter(){
        for (pos, qual) in methyl_info.iter(){
            aggregated_methyl_pos_dict.entry(*pos).or_insert(Vec::new()).push(*qual);
        }
    }
    let mut mod_score_dict: HashMap<usize, f32> = HashMap::new();
    for (pos, qual_list) in aggregated_methyl_pos_dict.iter(){
        let mut hyper_methyl = 0;
        let mut hypo_methyl = 0;
        for qual in qual_list.iter(){
            if *qual > min_methyl_threshold {
                hyper_methyl += 1;
            }else if *qual < (1.0 - min_methyl_threshold) {
                hypo_methyl += 1;
            }else{
                continue
            }
        }
 
        let mod_score = if hyper_methyl + hypo_methyl == 0 { 0.0 } else { hyper_methyl as f32 / (hyper_methyl + hypo_methyl) as f32 };
        mod_score_dict.insert(*pos, mod_score);
    }
    mod_score_dict
}


pub fn start (bam_records: HashMap<String, Record>, read_coordinates: &HashMap<String, (u64, u64)>) -> HashMap<String,HashMap<usize, f32>> {
    println!("Processing Methylation Signals from BAM file");
    // println!("Read coordinates: {:?}", read_coordinates);
    let mut read_name_set = HashSet::new();
    let mut methyl_all_reads: HashMap<String,HashMap<usize, f32>>= HashMap::new(); // extract directly from reference map
    let mut strand_dict:HashMap<String, bool> = HashMap::new();
    for (read_name, r) in bam_records.iter(){
        if r.is_unmapped(){
            continue
        }
        if r.is_supplementary(){
            continue
        }
        // let read_name = String::from_utf8_lossy(&r.qname()).to_string();
        if !read_coordinates.contains_key(read_name){
            println!("Read name not found in read coordinates: {:?}", read_name);
            continue
        }
        let (read_start, read_end) = read_coordinates.get(read_name).unwrap().clone();
        strand_dict.insert(read_name.clone(), r.is_reverse());
        read_name_set.insert(read_name.clone());

        let methyl_pos_dict = get_methylation_read(&r, read_start as usize, read_end as usize, 'm');
        methyl_all_reads.insert(read_name.clone(), methyl_pos_dict);
    }
    methyl_all_reads
}