version 1.0

workflow MHCAnnotate {
    input {
        File hap1_bam
        File hap1_bai
        File hap2_bam
        File hap2_bai
        String locus
        String prefix
    }

    call SubsetBam as hap1 {
        input:
            whole_genome_bam = hap1_bam,
            whole_genome_bai = hap1_bai,
            locus = locus,
            prefix = prefix + ".hap1"
    }

    call SubsetBam as hap2 {
        input:
            whole_genome_bam = hap2_bam,
            whole_genome_bai = hap2_bai,
            locus = locus,
            prefix = prefix + ".hap2"
    }

    call annotate as annot1 {
        input:
            fasta = hap1.output_fasta,
    }

    call annotate as annot2 {
        input:
            fasta = hap2.output_fasta,
    }
    output {
        Array[File] hap1_annotation = annot1.annotation
        Array[File] hap2_annotation = annot2.annotation
    }
}

task annotate {
    input {
        File fasta
    }

    command <<<
        set -euxo pipefail

        mkdir -p contigs output

        # Split multi-contig FASTA into one file per contig (mhca only uses the first record)
        awk '/^>/ {
            if (f) close(f)
            name = substr($1, 2)
            gsub(/[^A-Za-z0-9_.-]/, "_", name)
            f = "contigs/" name ".fa"
            print > f
            next
        }
        { print > f }' ~{fasta}

        for fa in contigs/*.fa; do
            mhca annotate "${fa}" output
        done
    >>>

    output {
        Array[File] annotation = glob("output/*.gff")
    }

    runtime {
        preemptible: 1
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/mhc_annotator:v1"
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
        prefix: "prefix for output bam, bai, and fasta file names"
    }

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{whole_genome_bam} ~{whole_genome_bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam
        samtools fasta ~{prefix}.bam > ~{prefix}.fa
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
        File output_fasta = "~{prefix}.fa"
    }
}
