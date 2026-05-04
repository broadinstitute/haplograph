version 1.0

workflow Extract_seq_from_asm_aligned_bam {
    input {
        Array[File] hap1_asm_bam
        Array[File] hap1_asm_bai
        Array[File] hap2_asm_bam
        Array[File] hap2_asm_bai
        Array[String] sample_id_list
        File reference_fa
        String prefix
        String locus
    }

    scatter (ind in range(length(hap1_asm_bam))){
        call get_truth_haplotypes {
            input:
                truth_hap1_bam = hap1_asm_bam[ind],
                truth_hap1_bai = hap1_asm_bai[ind],
                truth_hap2_bam = hap2_asm_bam[ind],
                truth_hap2_bai = hap2_asm_bai[ind],
                locus = locus,
                prefix = sample_id_list[ind]
        }
    }

    call MergeFastaFile {
        input:
            fastas = get_truth_haplotypes.fasta_file,
            prefix = prefix
    }

    output {
        File fastafile = MergeFastaFile.merged_fasta
    }
}


task get_truth_haplotypes {
    input {
        File truth_hap1_bam
        File truth_hap1_bai
        File truth_hap2_bam
        File truth_hap2_bai
        String locus
        String prefix
    }

    command <<<
        set -euxo pipefail

        /haplograph/target/release/haplograph extract -b ~{truth_hap1_bam} \
                                                -s ~{prefix} \
                                                -o ~{prefix}.hap1 \
                                                -l ~{locus} 

        /haplograph/target/release/haplograph extract -b ~{truth_hap2_bam} \
                                                -s ~{prefix} \
                                                -o ~{prefix}.hap2 \
                                                -l ~{locus}

        cat ~{prefix}.*.fasta > ~{prefix}.~{locus}.fasta
        
    >>>

    output {
        File fasta_file = "~{prefix}.~{locus}.fasta"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:v0.1.0"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }
}


struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

struct DataTypeParameters {
    Int num_shards
    String map_preset
}


task MergeFastaFile{
    input{
        Array[File] fastas
        String prefix
    }
    command <<<
        set -euxo pipefail

        # if fasta files are gzipped, unzip them first
        for fasta in ~{sep=" " fastas}; do
            if [[ "$fasta" == *.gz ]]; then
                gunzip -c "$fasta"
            else
                cat "$fasta"
            fi
        done > ~{prefix}.fasta

        gzip ~{prefix}.fasta
    >>>
    
    output{
        File merged_fasta="~{prefix}.fasta.gz"
        
    }

    Int disk_size = 100

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:v0.1.0"
    }
}


