version 1.0
workflow Locityper {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File reference_fa
        File reference_fai
        File alleles_fasta
        File? count_file
        String prefix
        String gcs_output_dir
        File gene_bed
    }
    if ! ~{defined(count_file)} {
        call Locityper_ref_construction {
            input:
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                prefix = prefix
        }

        call save_to_cloud {
            input:
                count_file = Locityper_ref_construction.output_count_file,
                gcs_output_dir = gcs_output_dir
        }

    }


    call SubsetBam {
        input:
            whole_genome_bam = whole_genome_bam,
            whole_genome_bai = whole_genome_bai,
            bed_file = gene_bed,
            prefix = prefix
    }

    call parseBed {
        input:
            bed = gene_bed,
            output_prefix = prefix
    }

    scatter (pair in zip(parseBed.locuslist, parseBed.genelist)) {
        String locus = pair.left
        String gene_name = pair.right
        call Locityper_genotyping {
            input:
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                gene_name = gene_name,
                locus = locus,
                alleles_fasta = alleles_fasta,
                align_bam = SubsetBam.output_bam,
                align_bai = SubsetBam.output_bai,
                count_file = select_first([count_file, Locityper_ref_construction.output_count_file]),
                prefix = prefix,

        }
    }

    output {

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

task Locityper_ref_construction {

    input {
        File reference_fa
        File reference_fai
        String prefix
    }


    command <<<
        set -eo pipefail

        jellyfish count --canonical --lower-count 2 --out-counter-len 2 --mer-len 25 \
            --threads 8 --size 3G --output HG38.counts.jf ~{reference_fa}
        

    >>>
    runtime {
        preemptible: 0
        memory: "8 GB"
        cpu: "2"
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/lr-locityper:latest"
    }
    output {
        File output_count_file = "HG38.counts.jf"
    }
}

task save_to_cloud {

    input {
        File count_file
        String gcs_output_dir
    }


    command <<<
        set -eo pipefail

        gsutil cp -r ~{count_file} ~{gcs_output_dir}
        
    >>>
    runtime {
        preemptible: 0
        memory: "8 GB"
        cpu: "2"
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/lr-hla_la:v1"
    }
    output {
        File output_count_file = "HG38.counts.jf"
    }
}



task Locityper_genotyping {

    input {
        File reference_fa
        File reference_fai
        String gene_name
        String locus
        File alleles_fasta
        File align_bam
        File align_bai
        File count_file
        String prefix
    }


    command <<<
        set -eo pipefail

        locityper target -d db -r ~{reference_fa} -j ~{count_file} \
            -l ~{gene_name} ~{locus} ~{alleles_fasta}

        locityper preproc -a ~{align_bam} -r ~{reference_fa} -j ~{count_file} -o bg/~{prefix}

        locityper genotype -a ~{align_bam} -d db -p bg/~{prefix} -o analysis/~{prefix}

        /usr/local/locityper/extra/into_csv.py \
            -i analysis/~{prefix}/* -o gts.csv

        /usr/local/locityper/extra/into_fasta.py -i gts.csv -d db -o out-dir


    >>>
    runtime {
        preemptible: 0
        memory: "8 GB"
        cpu: "2"
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/lr-locityper:latest"
    }
    output {
        File output_fasta = glob("out-dir/*")
    }
}

task parseBed {

    input {
        File bed
        String output_prefix

        Int? preemptible_tries
    }


    command <<<
        set -eo pipefail

        python - --bed_file ~{bed} \
                 --output_file ~{output_prefix} \
                 <<-'EOF'
        import gzip
        import argparse


        def import_bed(bed_file):
            locus_list = []
            gene_list = []
            if bed_file.endswith("gz"):
                with gzip.open(bed_file, "r") as f:
                    for line in f:
                        itemlist = line.strip().split("\t")
                        locus = "%s:%d-%d" % (itemlist[0], int(itemlist[1]), int(itemlist[2]))
                        locus_list.append(locus)
                        gene_list.append(itemlist[3])
            else:
                with open(bed_file, "r") as f:
                    for line in f:
                        itemlist = line.strip().split("\t")
                        locus = "%s:%d-%d" % (itemlist[0], int(itemlist[1]), int(itemlist[2]))
                        locus_list.append(locus)
                        gene_list.append(itemlist[3])
            return (locus_list, gene_list)


        def write_file(content, output_file):
            with open(output_file, "w") as f:
                for item in content:
                    f.write(f"{item}\n")

        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--bed_file',
                                type=str)

            parser.add_argument('--output_file',
                                type=str)

            args = parser.parse_args()

            locus_list, gene_list = import_bed(args.bed_file)
            write_file(locus_list, args.output_file + "_locus.txt")
            write_file(gene_list, args.output_file + "_gene.txt")


        if __name__ == "__main__":
            main()
        EOF

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/slee/kage-lite:pr_29"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }

    output {
        Array[String] locuslist = read_lines("~{output_prefix}_locus.txt")
        Array[String] genelist= read_lines("~{output_prefix}_gene.txt")
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