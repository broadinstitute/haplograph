version 1.0

workflow Haplograph_eval {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File reference_fa
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

        call haplograph {
            input:
                bam = SubsetBam.subsetbam,
                bai = SubsetBam.subsetbai,
                reference_fa = reference_fa,
                prefix = prefix + "_" + gene_name + "_" + desiredCoverage,
                locus = locus,
                windowsize = windowsize
        }






        
    }


    output {
        Array[Float] bam_coverage = CalculateCoverage.coverage
        Array[Array[File]] gfa = haplograph.graph_file
        Array[Array[File]] fasta = haplograph.asm_file
        Array[Array[File]] germline_vcf = haplograph.germline_vcf_file
        Array[Array[File]] somatic_vcf = haplograph.somatic_vcf_file
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


