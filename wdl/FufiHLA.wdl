version 1.0
workflow FuFiHLA {
    input {
        File whole_genome_bam
        File whole_genome_bai
        String prefix
        File gene_bed
    }

    call convertfasta {
        input:
            whole_genome_bam = whole_genome_bam,
            whole_genome_bai = whole_genome_bai,
            bed_file = gene_bed,
            prefix = prefix
    }

    call FuFiHLA_genotype {
        input:
            input_fa = convertfasta.output_fasta,
            prefix = prefix
    }

    output {
        FuFiHLAOutputs fufihla_output = FuFiHLA_genotype.outputs
    }
}

struct FuFiHLAOutputs {
    Array[File] HLA_A
    Array[File] HLA_B
    Array[File] HLA_C
    Array[File] HLA_DQA1
    Array[File] HLA_DQB1
    Array[File] HLA_DRB1
}

task convertfasta {

    input {
        File whole_genome_bam
        File whole_genome_bai
        File bed_file
        String prefix
    }

    parameter_meta {
        whole_genome_bam: {
            description: "bam to subset",
            localization_optional: true
        }
        whole_genome_bai:    "index for bam file"
        bed_file:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
    }


    command <<<
        set -eo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{whole_genome_bam} ~{whole_genome_bai} -L ~{bed_file} > ~{prefix}.bam
        samtools index ~{prefix}.bam
        samtools fasta ~{prefix}.bam > ~{prefix}.fa
        

    >>>
    runtime {
        preemptible: 0
        memory: "8 GB"
        cpu: "2"
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
    }
    output {
        File output_fasta = " ~{prefix}.fa"
    }
}

task FuFiHLA_genotype {
    input {
        File input_fa
        String prefix
        String extra_arg = "--hifi"
    }

    command <<<
        set -euxo pipefail

        fufihla --fa ~{input_fa} \
                --out ~{prefix} \
                ~{extra_arg}
        ls ~{prefix}
        ls ~{prefix}/consensus

    >>>

    output {
        FuFiHLAOutputs outputs = {
            "HLA_A": glob("~{prefix}/consensus/HLA-A*.fa"),
            "HLA_B": glob("~{prefix}/consensus/HLA-B*.fa"),
            "HLA_C": glob("~{prefix}/consensus/HLA-C*.fa"),
            "HLA_DQA1": glob("~{prefix}/consensus/HLA-DQA1*.fa"),
            "HLA_DQB1": glob("~{prefix}/consensus/HLA-DQB1*.fa"),
            "HLA_DRB1": glob("~{prefix}/consensus/HLA-DRB1*.fa"),
        }
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/lr-fufihla:latest"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 200 SSD"
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

