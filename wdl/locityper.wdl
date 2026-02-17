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
    if (! defined(count_file)) {
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
                read_fq = SubsetBam.output_fastq,
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
        samtools fastq ~{prefix}.bam > ~{prefix}.fq
        

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
        File output_fastq = "~{prefix}.fq"
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
        memory: "64 GB"
        cpu: "16"
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


task GenerateDB {
    input {
        File reference
        File reference_index
        File counts_jf
        String locus_name
        String locus_coordinates
        File alleles_fa
    }

    Int disk_size = 10
    String output_tar = locus_name + "db.tar.gz"

    command <<<
        set -euxo pipefail

        gunzip -c ~{reference} > reference.fa
        samtools faidx reference.fa

        locityper add -d ~{locus_name}.db \
            -r reference.fa \
            -j ~{counts_jf} \
            -l ~{locus_name} ~{locus_coordinates} ~{alleles_fa}

        find ~{locus_name}.db -type f -exec ls -lah {} \;
        
        echo "compressing DB"
        tar -czf ~{output_tar} ~{locus_name}.db
        echo "done compressing DB"
    >>>

    output {
        File db_tar = output_tar
    }

    runtime {
        memory: "8 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/lr-locityper:latest"
    }
}


task LocityperPreprocessAndGenotype {
    input {
        File input_fq
        File counts_file
        File reference
        File reference_index
        File db_targz
        String sample_id
        Array[String] locus_names

        Int locityper_n_cpu
        Int locityper_mem_gb

        String docker = "eichlerlab/locityper:0.19.1"
    }

    Int disk_size = 1 + 1*length(locus_names) + 2*ceil(size(select_all([input_fq, counts_file, reference, reference_index, db_targz]), "GiB"))
    String output_tar = sample_id + ".locityper.tar.gz"

    command <<<
        set -euxo pipefail
        
        df -h

        gunzip -c ~{reference} > reference.fa
        samtools faidx reference.fa

        nthreads=$(nproc)
        echo "using ${nthreads} threads"

        mkdir -p locityper_prepoc
        locityper preproc -i ~{input_fq} \
            -j ~{counts_file} \
            -@ ${nthreads} \
            --technology hifi \
            -r reference.fa \
            -o locityper_prepoc

        mkdir -p db
        tar --strip-components 1 -C db -xvzf ~{db_targz}

        mkdir -p out_dir
        locityper genotype -i ~{input_fq} \
            -d db \
            -p locityper_prepoc \
            -@ ${nthreads} \
            --debug 2 \
            --subset-loci ~{sep=" " locus_names} \
            -o out_dir

        tar -czf ~{output_tar} out_dir
    >>>

    runtime {
        memory: "~{locityper_mem_gb} GB"
        cpu: locityper_n_cpu
        disks: "local-disk ~{disk_size} HDD"
        preemptible: 0
        docker: docker
    }

    output {
        File genotype_tar = output_tar
    }
}


task Locityper_genotyping {

    input {
        File reference_fa
        File reference_fai
        String gene_name
        String locus
        File alleles_fasta
        File read_fq
        File count_file
        String prefix
        Int threads = 8
    }

    command <<<
        set -eo pipefail

        mkdir -p bg/~{prefix}

        locityper preproc -i ~{read_fq} \
                          -r ~{reference_fa} \
                          -@ ~{threads} \
                          -j ~{count_file} \
                          --technology hifi \
                          -o bg/~{prefix}

        mkdir -p analysis/~{prefix}

        locityper genotype -i ~{read_fq} \
                           -d db \
                           -p bg/~{prefix} \
                           -@ ~{threads} \
                           --subset-loci ~{gene_name} \
                           -o analysis/~{prefix}

        /usr/local/locityper/extra/into_csv.py \
            -i analysis/~{prefix}/* -o gts.csv

        mkdir -p out-dir
        /usr/local/locityper/extra/into_fasta.py -i gts.csv -d db -o out-dir


    >>>
    runtime {
        preemptible: 0
        memory: "64 GB"
        cpu: "16"
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