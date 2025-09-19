version 1.0

workflow Haplograph_eval {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File truth_asm_1_bam
        File truth_asm_1_bai
        File truth_asm_2_bam
        File truth_asm_2_bai
        File reference_fa
        String prefix
        String locus
        Int windowsize
        Array[Int] desiredCoverages
    }

    call CalculateCoverage {
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
            locus = locus,
            prefix = prefix
    }

    call get_truth_haplotypes {
        input:
            truth_hap1_bam = truth_asm_1_bam,
            truth_hap1_bai = truth_asm_1_bai,
            truth_hap2_bam = truth_asm_2_bam,
            truth_hap2_bai = truth_asm_2_bai,
            reference_fa = reference_fa,
            sampleid = prefix,
            locus = locus,
            prefix = prefix,
            extra_arg = "",
            windowsize = 1000,

    }

    scatter (desiredCoverage in desiredCoverages) {
        
        call downsampleBam {input:
            input_bam = CalculateCoverage.subsetbam,
            input_bam_bai = CalculateCoverage.subsetbai,
            basename = prefix,
            desiredCoverage = desiredCoverage,
            currentCoverage = CalculateCoverage.coverage,
            preemptible_tries = 0
        }

        call haplograph_asm {
            input:
                bam = downsampleBam.downsampled_bam,
                bai = downsampleBam.downsampled_bai,
                reference_fa = reference_fa,
                prefix = prefix,
                locus = locus,
                windowsize = windowsize
        }

        call haplograph_eval {
            input:
                truth_fasta = get_truth_haplotypes.fasta_file,
                query_fasta = haplograph_asm.asm_file,
                seq_number = 2,
                prefix = prefix + locus
        }
    }

    output {
        Float bam_coverage = CalculateCoverage.coverage
        Array[File] gfa = haplograph_asm.graph_file
        Array[File] fasta = haplograph_asm.asm_file
        Array[File] eval = haplograph_eval.qv_scores
    }
}

task haplograph_asm {
    input {
        File bam
        File bai
        File reference_fa
        String prefix
        String locus
        Int windowsize
        String extra_arg = ""
    }

    command <<<
        set -euxo pipefail
        /haplograph/target/release/haplograph haplotype -a ~{bam} \
                                                        -r ~{reference_fa} \
                                                        -s ~{prefix} \
                                                        -o ~{prefix} \
                                                        -l ~{locus} \
                                                        -w ~{windowsize} \
                                                        -d gfa \
                                                        ~{extra_arg}
        /haplograph/target/release/haplograph assemble -m -g ~{prefix}.gfa -o ~{prefix}.fasta
        
    >>>

    output {
        File graph_file = "~{prefix}.gfa"
        File asm_file = "~{prefix}.fasta"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograp:v1"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }
}


task haplograph_eval {
    input {
        File truth_fasta
        File query_fasta
        Int seq_number
        String prefix
    }

    command <<<
        set -euxo pipefail

        /haplograph/target/release/haplograph evaluate -t ~{truth_fasta} \
                                                        -q ~{query_fasta} \
                                                        -s ~{seq_number} \
                                                        -o ~{prefix}.tsv
        
    >>>

    output {
        File qv_scores = "~{prefix}.tsv"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograp:v1"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }
}

task get_truth_haplotypes_from_annotation {
    input {
        File truth_asm_1
        File truth_asm_annotation_1
        File truth_asm_2
        File truth_asm_annotation_2
        File reference_fa
        String gene_name
        String prefix
    }

    command <<<
        set -euxo pipefail

        python - --truth_fasta_1 ~{truth_asm_1} \
                 --truth_annotation_1 ~{truth_asm_annotation_1} \
                 --truth_fasta_2 ~{truth_asm_2} \
                 --truth_annotation_2 ~{truth_asm_annotation_2} \
                 --gene_name ~{gene_name} \
                 --output_filename ~{prefix}.~{gene_name}.truth.fasta \
                 --output_header ~{gene_name} \
                 <<-'EOF'
        import pysam
        import gzip
        import argparse

        def get_hla_region(annotation_file, name):
            intervals = []
            with gzip.open(annotation_file, "rt") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    itemlist = line.strip().split("\t")
                    if itemlist[2] == "gene":
                        contig = itemlist[0]
                        start = int(itemlist[3])
                        end = int(itemlist[4])
                        for item in itemlist[8].split(";"):
                            if "gene_name" in item:
                                gene_name = item.split(" ")[2]
                                if name in gene_name:
                                    intervals.append((contig, start, end))
            return intervals  

        def get_truth_seq(fastafile, annotation_file, gene_name):
            fasta_records = pysam.Fastafile(fastafile)
            chromosomes = fasta_records.references
            truth_seq = []
            intervals = get_hla_region(annotation_file, gene_name)
            for chromosome in chromosomes:
                seq = fasta_records.fetch(chromosome)
                for (contig_name, start, end) in intervals:
                    if contig_name == chromosome:
                        truth_seq.append(seq[start:end])
            return truth_seq  


        def write_fasta(fasta_file, header, truth_seq):
            with open(fasta_file, "w") as f:
                for i, t_seq in enumerate(truth_seq):
                    f.write(">%s_Hap_%d\n" % (header, i+1))
                    char_length = 60
                    line_number = len(t_seq) // char_length

                    for i in range(line_number):
                        f.write(t_seq[i*char_length:(i+1)*char_length] + "\n")
                    f.write(t_seq[line_number*char_length:] + "\n")



        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--truth_fasta_1',
                                type=str)

            parser.add_argument('--truth_annotation_1',
                                type=str)

            parser.add_argument('--truth_fasta_2',
                                type=str)

            parser.add_argument('--truth_annotation_2',
                                type=str)

            parser.add_argument('--gene_name',
                                type=str)

            parser.add_argument('--output_filename',
                                type=str)
            parser.add_argument('--output_header',
                            type=str)

            args = parser.parse_args()

            truth_seq_1 = get_truth_seq(args.truth_fasta_1, args.truth_annotation_1, args.gene_name)
            truth_seq_2 = get_truth_seq(args.truth_fasta_2, args.truth_annotation_2, args.gene_name)
            truth_seq = truth_seq_1 + truth_seq_2
            write_fasta(args.output_filename, args.output_header, truth_seq)

        if __name__ == "__main__":
            main()
        EOF
        
    >>>

    output {
        File fasta_file = "~{prefix}.~{gene_name}.truth.fasta"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/slee/kage-lite:pr_29"
        memory: "8 GB"
        cpu: 2
        disks: "local-disk 100 SSD"
    }
}


task get_truth_haplotypes {
    input {
        File truth_hap1_bam
        File truth_hap1_bai
        File truth_hap2_bam
        File truth_hap2_bai
        File reference_fa
        String sampleid
        String locus
        String prefix
        String extra_arg
        Int windowsize
    }

    command <<<
        set -euxo pipefail

        /haplograph/target/release/haplograph haplotype -a ~{truth_hap1_bam} \
                                                -r ~{reference_fa} \
                                                -s ~{sampleid} \
                                                -o ~{prefix}.truth1 \
                                                -l ~{locus} \
                                                -w ~{windowsize} \
                                                -m 0 \
                                                -f 0 \
                                                -d gfa \
                                                ~{extra_arg}

        /haplograph/target/release/haplograph assemble -g ~{prefix}.truth1.gfa -o ~{prefix}.truth1.fasta

        /haplograph/target/release/haplograph haplotype -a ~{truth_hap2_bam} \
                                                -r ~{reference_fa} \
                                                -s ~{sampleid} \
                                                -o ~{prefix}.truth2 \
                                                -l ~{locus} \
                                                -w ~{windowsize} \
                                                -m 0 \
                                                -f 0 \
                                                -d gfa \
                                                ~{extra_arg}
        /haplograph/target/release/haplograph assemble -g ~{prefix}.truth2.gfa -o ~{prefix}.truth2.fasta
        cat ~{prefix}.truth1.fasta ~{prefix}.truth2.fasta > ~{prefix}.truth.fasta
        
    >>>

    output {
        File fasta_file = "~{prefix}.truth.fasta"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograp:v1"
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

task CalculateCoverage {

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
        String prefix = "subset"

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
        mem_gb:             10,
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


task downsampleBam {

    input {
        File input_bam
        File input_bam_bai
        String basename
        Int desiredCoverage
        Float currentCoverage
        Int? preemptible_tries
    }

    meta {
        description: "Uses Picard to downsample to desired coverage based on provided estimate of coverage."
    }

    parameter_meta {
    }

    Float scalingFactor = desiredCoverage / currentCoverage


    command <<<
        set -eo pipefail

        gatk DownsampleSam -I ~{input_bam} -O ~{basename}_~{desiredCoverage}x.bam -R 7 -P ~{scalingFactor} -S ConstantMemory --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true


    >>>
    runtime {
        preemptible: select_first([preemptible_tries, 5])
        memory: "8 GB"
        cpu: "2"
        disks: "local-disk 500 HDD"
        docker: "us.gcr.io/broad-gatk/gatk"
    }
    output {
        File downsampled_bam = "~{basename}_~{desiredCoverage}x.bam"
        File downsampled_bai = "~{basename}_~{desiredCoverage}x.bai"
    }
}
