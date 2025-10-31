version 1.0

workflow Haplograph_eval {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File truth_asm_1
        File truth_asm_1_annotation
        File truth_asm_2
        File truth_asm_2_annotation
        File truth_vcf
        File? truth_vcf_tbi
        File reference_fa
        String prefix
        File gene_bed
        Int windowsize
        Array[Int] desiredCoverages
        Int hifiasm_mem
        Int hifiasm_thread
    }



    call parseBed {
        input:
            bed = gene_bed,
            output_prefix = prefix
    }

    scatter (pair in zip(parseBed.locuslist, parseBed.genelist)) {
        String locus = pair.left
        String gene_name = pair.right

        call CalculateCoverage {
            input:
                bam = whole_genome_bam,
                bai = whole_genome_bai,
                locus = locus,
                prefix = prefix + "_" + gene_name
        }

        call get_truth_haplotypes_from_annotation {
            input:
                truth_asm_1 = truth_asm_1,
                truth_asm_annotation_1 = truth_asm_1_annotation,
                truth_asm_2 = truth_asm_2,
                truth_asm_annotation_2 = truth_asm_2_annotation,
                reference_fa = reference_fa,
                gene_name = gene_name,
                prefix = prefix + "_" + gene_name
        }



        scatter (desiredCoverage in desiredCoverages) {
            
            call downsampleBam {input:
                input_bam = CalculateCoverage.subsetbam,
                input_bam_bai = CalculateCoverage.subsetbai,
                basename = prefix + "_" + gene_name,
                desiredCoverage = desiredCoverage,
                currentCoverage = CalculateCoverage.coverage,
                preemptible_tries = 0
            }

            call haplograph {
                input:
                    bam = downsampleBam.downsampled_bam,
                    bai = downsampleBam.downsampled_bai,
                    reference_fa = reference_fa,
                    prefix = prefix + "_" + gene_name + "_" + desiredCoverage,
                    locus = locus,
                    windowsize = windowsize
            }

            call haplograph_eval {
                input:
                    truth_fasta = get_truth_haplotypes_from_annotation.fasta_file,
                    query_fasta = haplograph.asm_file,
                    prefix = prefix + "_" + gene_name + "_" + desiredCoverage,
            }

            call hifiasm_asm {
                input:
                    bam = downsampleBam.downsampled_bam,
                    bai = downsampleBam.downsampled_bai,
                    prefix = prefix + "_" + gene_name + "_" + desiredCoverage,
                    num_cpus = hifiasm_thread,
                    mem_gb = hifiasm_mem
            }

            call haplograph_eval as hifiasm_eval {
                input:
                    truth_fasta = get_truth_haplotypes_from_annotation.fasta_file,
                    query_fasta = hifiasm_asm.asm_file,
                    prefix = prefix + "_" + gene_name + "_" + desiredCoverage,
            }

            call Vcfdist as VCFdist_germline {
                input:
                    sample = prefix,
                    eval_vcf = haplograph.vcf_file,
                    truth_vcf = truth_vcf,
                    locus = locus,
                    reference_fasta = reference_fa,
                    coverage = desiredCoverage,
                    genename = gene_name,
                    extra_args = ""
            }
        }
    }


    output {
        Array[Float] bam_coverage = CalculateCoverage.coverage
        Array[Array[File]] gfa = haplograph.graph_file
        Array[Array[File]] fasta = haplograph.asm_file
        Array[Array[File]] vcf = haplograph.vcf_file
        Array[Array[File]] haplograph_eval_result = haplograph_eval.qv_scores
        Array[Array[File]] hifiasm_eval_result = hifiasm_eval.qv_scores
        Array[Array[VcfdistOutputs]] vcfdist_summary = VCFdist_germline.outputs
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
        String extra_arg = ""
    }

    command <<<
        set -euxo pipefail
        /haplograph/target/release/haplograph haplograph -a ~{bam} \
                                                        -r ~{reference_fa} \
                                                        -s ~{prefix} \
                                                        -o ~{prefix} \
                                                        -l ~{locus} \
                                                        -w ~{windowsize} \
                                                        -m ~{minimal_supported_reads} \
                                                        -d gfa \
                                                        ~{extra_arg}
        
        
    >>>

    output {
        File graph_file = "~{prefix}.gfa"
        File asm_file = "~{prefix}.fasta"
        File vcf_file = "~{prefix}.vcf.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:v2"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}


task haplograph_eval {
    input {
        File truth_fasta
        File query_fasta
        String prefix
    }

    command <<<
        set -euxo pipefail
        /haplograph/target/release/haplograph evaluate -t ~{truth_fasta} \
                                                        -q ~{query_fasta} \
                                                        -s 2 \
                                                        -o ~{prefix}.tsv
        
    >>>

    output {
        File qv_scores = "~{prefix}.tsv"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:v2"
        memory: "16 GB"
        cpu: 4
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
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/haplograph:v2"
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
        if awk "BEGIN{ if((~{scalingFactor}) < 1.0) exit 0; else exit 1 }"; then
            gatk DownsampleSam -I ~{input_bam} -O ~{basename}_~{desiredCoverage}x.bam -R 7 -P ~{scalingFactor} -S ConstantMemory --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true
        else
            mv ~{input_bam} ~{basename}_~{desiredCoverage}x.bam
            mv ~{input_bam_bai} ~{basename}_~{desiredCoverage}x.bai
            echo "Total Coverage is lower than desiredCoverage"
        fi

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

struct VcfdistOutputs {
    File summary_vcf
    File precision_recall_summary_tsv
    File precision_recall_tsv
    File query_tsv
    File truth_tsv
    File phasing_summary_tsv
    File switchflips_tsv
    File superclusters_tsv
    File phase_blocks_tsv
}

task Vcfdist {
    input {
        String truth_sample
        String sample
        String genename
        String coverage
        File eval_vcf
        File? eval_vcf_index
        File truth_vcf
        File? truth_vcf_index
        String locus
        File reference_fasta
        String? extra_args
        Int verbosity = 1

        Int disk_size_gb = ceil(size(truth_vcf, "GiB") + 10)
        Int mem_gb = 16
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail
        bcftools index -t ~{truth_vcf}
        bcftools view -s ~{truth_sample} -r ~{locus} ~{truth_vcf} -Oz -o ~{sample}.~{locus}.base.vcf.gz
        bcftools index -t ~{sample}.~{locus}.base.vcf.gz

        bcftools index -t ~{eval_vcf}
        bcftools filter -e 'INFO/SOMATIC=1' ~{eval_vcf} -Oz -o ~{sample}.filtered.query.vcf.gz
        bcftools index -t ~{sample}.filtered.query.vcf.gz

        vcfdist \
            ~{sample}.filtered.query.vcf.gz \
            ~{sample}.~{locus}.base.vcf.gz \
            ~{reference_fasta} \
            -v ~{verbosity} \
            ~{extra_args}

        for tsv in $(ls *.tsv); do mv $tsv ~{sample}.~{genename}.~{coverage}.$tsv; done
        mv summary.vcf ~{sample}.~{genename}.~{coverage}.summary.vcf
    >>>

    output {
        VcfdistOutputs outputs = {
            "summary_vcf": "~{sample}.~{genename}.~{coverage}.summary.vcf",
            "precision_recall_summary_tsv": "~{sample}.~{genename}.~{coverage}.precision-recall-summary.tsv",
            "precision_recall_tsv": "~{sample}.~{genename}.~{coverage}.precision-recall.tsv",
            "query_tsv": "~{sample}.~{genename}.~{coverage}.query.tsv",
            "truth_tsv": "~{sample}.~{genename}.~{coverage}.truth.tsv",
            "phasing_summary_tsv": "~{sample}.~{genename}.~{coverage}.phasing-summary.tsv",
            "switchflips_tsv": "~{sample}.~{genename}.~{coverage}.switchflips.tsv",
            "superclusters_tsv": "~{sample}.~{genename}.~{coverage}.superclusters.tsv",
            "phase_blocks_tsv": "~{sample}.~{genename}.~{coverage}.phase-blocks.tsv"
        }
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/vcfdist:v1"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}

task hifiasm_asm{
    input{
        File bam
        File bai
        String prefix
        Int num_cpus
        Int mem_gb
    }

    Int disk_size = 10 + ceil(2 * size(bam, "GiB"))

    command <<<

        set -euxo pipefail

        samtools fastq ~{bam} > ~{prefix}.fastq

        /truvari/hifiasm-0.25.0/hifiasm -o ~{prefix} -t 4 ~{prefix}.fastq
        awk '/^S/{print ">"$2;print $3}' ~{prefix}.bp.hap1.p_ctg.gfa > ~{prefix}.bp.hap1.p_ctg.fa
        awk '/^S/{print ">"$2;print $3}' ~{prefix}.bp.hap2.p_ctg.gfa > ~{prefix}.bp.hap2.p_ctg.fa

        cat ~{prefix}.bp.hap1.p_ctg.fa ~{prefix}.bp.hap2.p_ctg.fa > ~{prefix}.hifiasm.fa
    >>>

    output{
        File assembly_hap1="~{prefix}.bp.hap1.p_ctg.fa"
        File assembly_hap2="~{prefix}.bp.hap2.p_ctg.fa"
        File asm_file = "~{prefix}.hifiasm.fa"
    }
    runtime {
        cpu: num_cpus
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/hifiasm:0.25.0"
    }    
}

