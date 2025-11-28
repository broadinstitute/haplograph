version 1.0

workflow HLA_LA {
    input {
        File whole_genome_bam
        File whole_genome_bai
        String gcs_input_dir
        File gene_bed
        String prefix
    }

    call SubsetBam {
        input:
            whole_genome_bam = whole_genome_bam,
            whole_genome_bai = whole_genome_bai,
            bed_file = gene_bed,
            prefix = prefix
    }

    call HLA_LA_typing {
        input:
            bam = SubsetBam.output_bam,
            bai = SubsetBam.output_bai,
            sample_id = prefix,
            gcs_input_dir = gcs_input_dir
    }



    output {
        Array[File] output_files = HLA_LA_typing.output_file
    }
}

task HLA_LA_build_graph {
    input {
        String gcs_output_dir
    }

    command <<<
        set -euxo pipefail

        mkdir graphs
        cd graphs
        wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz
        tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz
        
        HLA-LA --action prepareGraph --PRG_graph_dir ./PRG_MHC_GRCh38_withIMGT

        gsutil rsync -r . ~{gcs_output_dir}


    >>>

    output {
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/lr-hla_la:v1"
        memory: "64 GB"
        cpu: 16
        disks: "local-disk 200 SSD"
    }
}

task HLA_LA_typing {
    input {
        File bam
        File bai
        String sample_id
        String gcs_input_dir
        String extra_arg = ""
    }

    command <<<
        set -euxo pipefail

        cd /usr/local/bin/HLA-LA/graphs
        gsutil rsync -r ~{gcs_input_dir} .

        cd /mnt/disks/cromwell_root/ 

        /usr/local/bin/HLA-LA/src/HLA-LA.pl --BAM ~{bam} \
                                            --graph PRG_MHC_GRCh38_withIMGT \
                                            --sampleID ~{sample_id} \
                                            --maxThreads 7 \
                                            --bwa_bin /usr/bin \
                                            --longReads pacbio \
                                            --workingDir .

    >>>

    output {
        Array[File] output_file = glob("*")
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/lr-hla_la:v1"
        memory: "64 GB"
        cpu: 16
        disks: "local-disk 200 SSD"
    }
}

task SubsetBam {

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
        File output_bam = "~{prefix}.bam"
        File output_bai = "~{prefix}.bam.bai"
        File output_fasta = "~{prefix}.fa"
    }
}
