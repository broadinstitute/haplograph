version 1.0

workflow StarPhase_HLA {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File database
        File ref_fa
        String locus
        String prefix
    }

    call SubsetBam {
        input:
            whole_genome_bam = whole_genome_bam,
            whole_genome_bai = whole_genome_bai,
            locus = locus,
            prefix = prefix
    }

    call StarPhase {
        input:
            database = database,
            reference_fasta = ref_fa,
            bam = SubsetBam.output_bam,
            bai = SubsetBam.output_bai,
            prefix = prefix
    }

    output {
        File output_files = StarPhase.starphase_output
    }
}

task StarPhase {
    input {
        File database
        File reference_fasta
        File bam
        File bai
        String prefix
    }

    String bam_basename = basename(bam)

    command <<<
        set -euxo pipefail

        # Cromwell localizes bam/bai under different paths; co-locate for tools that expect adjacent index
        ln -s ~{bam} ~{bam_basename}
        ln -s ~{bai} ~{bam_basename}.bai

        pbstarphase diplotype \
            --database ~{database} \
            --reference ~{reference_fasta} \
            --bam ~{bam_basename} \
            --output-calls ~{prefix}.json
    >>>

    output {
        File starphase_output = "~{prefix}.json"
    }

    runtime {
        preemptible: 1
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/starphase:2.0.1"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 200 SSD"
    }
}

task SubsetBam {

    input {
        File whole_genome_bam
        File whole_genome_bai
        String locus
        String prefix
    }

    parameter_meta {
        whole_genome_bam: {
            description: "bam to subset",
            localization_optional: true
        }
        whole_genome_bai: "index for bam file"
        locus: "genomic locus to select (e.g. chr6:28510120-33480577)"
        prefix: "prefix for output bam and bai file names"
    }

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{whole_genome_bam} ~{whole_genome_bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam
    >>>

    runtime {
        preemptible: 1
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
    }

    output {
        File output_bam = "~{prefix}.bam"
        File output_bai = "~{prefix}.bam.bai"
    }
}
