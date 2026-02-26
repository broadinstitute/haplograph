version 1.0

workflow Haplopan_eval {
    input {
        File whole_genome_bam
        File whole_genome_bai
        Array[File] reference_asm_hap1_bam
        Array[File] reference_asm_hap1_bai
        Array[File] reference_asm_hap2_bam
        Array[File] reference_asm_hap2_bai
        String prefix
        File gene_bed
        Int windowsize
    }

    call parseBed {
        input:
            bed = gene_bed,
            output_prefix = prefix
    }

    scatter (pair in zip(parseBed.locuslist, parseBed.genelist)) {
        String locus = pair.left
        String gene_name = pair.right

        call SubsetBam {
            input:
                bam = whole_genome_bam,
                bai = whole_genome_bai,
                locus = locus,
                prefix = prefix + "_" + gene_name
        }

        call BuildPangenome{
            input:
                hap1_aligned_bam = reference_asm_hap1_bam,
                hap1_aligned_bai = reference_asm_hap1_bai,
                hap2_aligned_bam = reference_asm_hap2_bam,
                hap2_aligned_bai = reference_asm_hap2_bai,
                locus = locus,
                prefix = prefix + "_" + gene_name,
                sample_id = prefix
        }

        call haplopan {
            input:
                bam = SubsetBam.subsetbam,
                bai = SubsetBam.subsetbai,
                pangenome_fa = BuildPangenome.extracted_fasta,
                prefix = prefix + "_" + gene_name,
        }

        
    }


    output {
        Array[Array[File]] gfa = haplopan.graph_file
        Array[File] fasta = haplopan.asm_file
    }
}

task haplopan {
    input {
        File bam
        File bai
        File pangenome_fa
        String prefix
        String sample_id
        String extra_arg = ""
    }

    command <<<
        set -euxo pipefail
        /haplograph/target/release/haplograph haplopan -p ~{pangenome_fa} \
                                                        -a ~{bam} \
                                                        -o ~{prefix} \
                                                        -s ~{sample_id}
                                                        ~{extra_arg}
        
        ls -l .
    >>>

    output {
        Array[File] graph_file = glob("*.gfa")
        File asm_file = "~{prefix}.final.fasta "
        Array[File] methyl_bed = glob("*.bed")
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:v4"
        memory: "32 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}


task haplograph {
    input {
        File bam
        File bai
        File reference_fa
        String prefix
        String locus
        Int windowsize
        Int minimal_supported_reads
        Int fold_threshold
        Float min_freq
        String extra_arg = ""
    }

    command <<<
        set -euxo pipefail
        /haplograph/target/release/haplograph haplograph -a ~{bam} \
                                                        -r ~{reference_fa} \
                                                        -s ~{prefix} \
                                                        -o ~{prefix} \
                                                        -l ~{locus} \
                                                        -v ~{min_freq} \
                                                        -m ~{minimal_supported_reads} \
                                                        -w ~{windowsize} \
                                                        -f gfa \
                                                        -c ~{fold_threshold}
                                                        ~{extra_arg}
        
        ls -l .
    >>>

    output {
        File graph_file = "~{prefix}.gfa"
        File asm_file = "~{prefix}.fasta"
        File germline_vcf_file = "~{prefix}.germline.vcf.gz"
        File somatic_vcf_file = "~{prefix}.somatic.vcf.gz"
        Array[File] methyl_bed = glob("*.bed")
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:v3"
        memory: "16 GB"
        cpu: 4
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

task SubsetBam {

    meta {
        description : "Subset a BAM file to a specified locus."
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        File bai
        String locus
        String prefix

        RuntimeAttr? runtime_attr_override
    }



    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam
        samtools depth -r ~{locus} ~{prefix}.bam | awk '{sum+=$3} END {print sum/NR}' > coverage.txt

    >>>

    output {
        Float coverage = read_float("coverage.txt")
        File subsetbam =  "~{prefix}.bam"
        File subsetbai = " ~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task BuildPangenome {

    meta {
        description : "Build a Pangenome Fasta file to a specified locus."
    }

    parameter_meta {
        hap1_aligned_bam: {
            description: "bam to subset",
            localization_optional: true
        }
        hap1_aligned_bai:    "index for bam file"

        hap2_aligned_bam: {
            description: "bam to subset",
            localization_optional: true
        }
        hap2_aligned_bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
        sample_id: "sample ID"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        Array[File] hap1_aligned_bam
        Array[File] hap1_aligned_bai
        Array[File] hap2_aligned_bam
        Array[File] hap2_aligned_bai
        String locus
        String prefix
        String sample_id

        RuntimeAttr? runtime_attr_override
    }

    Array[File] aligned_bams = flatten([hap1_aligned_bam, hap2_aligned_bam])

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        index=0

        for bam in ~{sep=" " aligned_bams}; do
            /haplograph/target/release/haplograph extract -b "$bam" -l "~{locus}" -o "${index}.fasta" -s "~{sample_id}"
            index=$((index + 1))
        done

        cat *.fasta > ~{prefix}.fasta

    >>>

    output {
        File extracted_fasta = " ~{prefix}.fasta"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:v4"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
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


