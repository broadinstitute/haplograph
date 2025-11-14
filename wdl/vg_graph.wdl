
version 1.0

workflow VG_extract_bam_haplograph {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File gfa_file
        File truth_asm_1
        File truth_asm_1_annotation
        File truth_asm_2
        File truth_asm_2_annotation
        File truth_vcf
        File gene_bed
        Int windowsize
        Array[Int] desiredCoverages
        Int hifiasm_mem
        Int hifiasm_thread
        File reference_fa

        String prefix
        String mhc_hg38_locus
        Int mhc_start_pos
        String mhc_path_name

        Int giraffe_threads
        String sample_id
        String giraffe_extra_arg
        String graph_ref_path
    }

    call SubsetBam {
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
            locus = mhc_hg38_locus,
            prefix = prefix
    }

    call vg_graph_index {
        input:
            graph = gfa_file,
            prefix = prefix
    }

    call vg_graph_alignment {
        input:
            graph_index = vg_graph_index.graph_index,
            graph_dist = vg_graph_index.graph_dist,
            graph_min = vg_graph_index.graph_min,
            graph_zipcode = vg_graph_index.graph_zipcode,
            input_fq = SubsetBam.subsetfq,
            thread_num = giraffe_threads,
            in_sample_name = sample_id,
            extra_args = giraffe_extra_arg,
            out_prefix = prefix + "_vg_align",
    }

    call vg_graph_surject {
        input:
            graph_index = vg_graph_index.graph_index,
            input_gam = vg_graph_alignment.gam_file,
            thread_num = giraffe_threads,
            in_sample_name = sample_id,
            extra_args = giraffe_extra_arg,
            out_prefix = prefix + "_vg_align",
            ref_path = graph_ref_path
    }
    
    call MHC_parseBed {
        input:
            bed = gene_bed,
            spos = mhc_start_pos,
            mhc_path_name = mhc_path_name,
            output_prefix = prefix
    }

    scatter (i in range(length(MHC_parseBed.locuslist))) {
            String locus = MHC_parseBed.locuslist[i]
            String hg38_locus = MHC_parseBed.reflocuslist[i]
            String gene_name = MHC_parseBed.genelist[i]

            call CalculateCoverage {
                input:
                    bam = vg_graph_surject.bam_file,
                    bai = vg_graph_surject.bam_index_file,
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
                        reference_fa = vg_graph_index.graph_reference_fa,
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

                call liftover_vcf {
                    input:
                        vcf = haplograph.germline_vcf_file,
                        current_contig = mhc_path_name,
                        target_contig = "chr6",
                        path_spos = 28510125,
                        prefix = prefix + "_" + gene_name + "_" + desiredCoverage + "liftover",
                }

                call Vcfdist as VCFdist_germline {
                    input:
                        sample = prefix,
                        eval_vcf = liftover_vcf.lifted_vcf,
                        truth_vcf = truth_vcf,
                        locus = hg38_locus,
                        reference_fasta = reference_fa,
                        coverage = desiredCoverage,
                        genename = gene_name
                }


            }
    }
    
    output {
        File graph_bam = vg_graph_surject.bam_file
        File graph_bai = vg_graph_surject.bam_index_file
        Array[Array[File]] haplograph_eval_result = haplograph_eval.qv_scores
        Array[Array[File]] hifiasm_eval_result = hifiasm_eval.qv_scores
        Array[Array[File]] vcfdist_call_summary = VCFdist_germline.precision_recall_summary_tsv
        Array[Array[File]] vcfdist_phasing_summary = VCFdist_germline.phasing_summary_tsv

    }
}

task SubsetBam {

    input {
        File bam
        File bai
        String locus
        String prefix
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
    }


    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam

        samtools fastq -F 0x900 -OT RG,BC,Mm,Ml -n ~{prefix}.bam > ~{prefix}.fastq
 
    >>>

    output {
        File subsetbam =  "~{prefix}.bam"
        File subsetbai = "~{prefix}.bam.bai"
        File subsetfq = "~{prefix}.fastq"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/vg:1.69.0"
        cpu: 1
        memory: "4GB"
        disks: "local-disk " + 100 + " HDD"
        preemptible: 0
    }

}



task vg_graph_index{
    input{
        File graph
        String prefix
    }

    command <<<

    set -euxo pipefail
    /vg/vg autoindex --workflow lr-giraffe \
                     --prefix ~{prefix} \
                     --gfa ~{graph}

    /vg/vg paths -x ~{graph} -R -F > ~{prefix}.ref.fa

    >>>

    output {
        File graph_index= "~{prefix}.giraffe.gbz"
        File graph_dist = "~{prefix}.dist"
        File graph_min = "~{prefix}.longread.withzip.min"
        File graph_zipcode = "~{prefix}.longread.zipcodes"
        File graph_reference_fa = "~{prefix}.ref.fa"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/vg:1.69.0"
        cpu: 1
        memory: "4GB"
        disks: "local-disk " + 100 + " HDD"
        preemptible: 0
    }
}

task vg_graph_alignment {

    input {
        File graph_index
        File graph_dist
        File graph_min
        File graph_zipcode
        File input_fq
        Int thread_num
        String read_type = "hifi"
        String in_sample_name
        String extra_args
        String out_prefix
    }
    Int half_cores = thread_num / 2

    command <<<
        set -euxo pipefail

        /vg/vg giraffe --parameter-preset ~{read_type} \
                       --progress --track-provenance \
                       -t ~{thread_num} \
                       -Z ~{graph_index} \
                       -d ~{graph_dist} \
                       -m ~{graph_min} \
                       -z ~{graph_zipcode} \
                       -f ~{input_fq} \
                       -o gam > ~{out_prefix}.gam     
        
    >>>

    output {
        File gam_file = "~{out_prefix}.gam"
    }

    #########################

    runtime {
        docker:  "us.gcr.io/broad-dsp-lrma/hangsuunc/vg:1.69.0"
        cpu: 16
        memory: "64GB"
        disks: "local-disk " + 200 + " HDD"
        preemptible: 0
    }
    
}

task vg_graph_surject {

    input {
        File graph_index
        Int thread_num
        File input_gam
        Boolean input_is_gam = true
        String read_type = "hifi"
        String in_sample_name
        String extra_args
        String out_prefix
        String ref_path
    }
    Int half_cores = thread_num / 2

    command <<<
        set -euxo pipefail

        /vg/vg surject \
                -n ~{ref_path} \
                -x ~{graph_index} \
                -t ~{thread_num} \
                --multimap \
                --bam-output ~{true="" false="--gaf-input" input_is_gam} \
                --sample ~{in_sample_name} \
                --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:hifi PU:unit1" \
                ~{extra_args} \
                ~{input_gam}| samtools sort --threads ~{half_cores} \
                                            -O BAM > ~{out_prefix}.bam
        samtools index -b ~{out_prefix}.bam
        
    >>>

    output {
        File bam_file = "~{out_prefix}.bam"
        File bam_index_file = "~{out_prefix}.bam.bai"
    }

    #########################

    runtime {
        docker:  "us.gcr.io/broad-dsp-lrma/hangsuunc/vg:1.69.0"
        cpu: 16
        memory: "64GB"
        disks: "local-disk " + 200 + " HDD"
        preemptible: 0
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

task Vcfdist {
    input {
        String truth_sample
        String sample
        String genename
        String coverage
        File eval_vcf
        File truth_vcf
        String locus
        File reference_fasta
        String extra_args = "--realign-truth --realign-query --verbosity 2 --largest-variant 1000"
        Int verbosity = 1

        Int disk_size_gb = ceil(size(truth_vcf, "GiB") + 10)
        Int mem_gb = 64
        Int cpu = 16
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        bcftools index -t ~{truth_vcf}
        bcftools view -s ~{truth_sample} -r ~{locus} ~{truth_vcf} -Oz -o ~{sample}.~{locus}.base.vcf.gz
        bcftools index -t ~{sample}.~{locus}.base.vcf.gz

        bcftools index -t ~{eval_vcf}
        bcftools filter -e 'INFO/SOMATIC=1' ~{eval_vcf} -Oz -o ~{sample}.filtered.query.vcf.gz
        bcftools +fixploidy ~{sample}.filtered.query.vcf.gz -- -f 2 > ~{sample}.filtered.query.fixploidy.vcf
        bcftools view ~{sample}.filtered.query.fixploidy.vcf -Oz -o ~{sample}.filtered.query.fixploidy.vcf.gz
        bcftools index -t ~{sample}.filtered.query.fixploidy.vcf.gz

        vcfdist \
            ~{sample}.filtered.query.fixploidy.vcf.gz \
            ~{sample}.~{locus}.base.vcf.gz \
            ~{reference_fasta} \
            -v ~{verbosity} \
            ~{extra_args}

        for tsv in $(ls *.tsv); do mv $tsv ~{sample}.~{genename}.~{coverage}.$tsv; done
        mv summary.vcf ~{sample}.~{genename}.~{coverage}.summary.vcf
    >>>

    output {
        File summary_vcf = "~{sample}.~{genename}.~{coverage}.summary.vcf"
        File precision_recall_summary_tsv = "~{sample}.~{genename}.~{coverage}.precision-recall-summary.tsv"
        File query_tsv = "~{sample}.~{genename}.~{coverage}.query.tsv"
        File truth_tsv = "~{sample}.~{genename}.~{coverage}.truth.tsv"
        File phasing_summary_tsv = "~{sample}.~{genename}.~{coverage}.phasing-summary.tsv"
        File switchflips_tsv = "~{sample}.~{genename}.~{coverage}.switchflips.tsv"
        File phase_blocks_tsv = "~{sample}.~{genename}.~{coverage}.phase-blocks.tsv"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/vcfdist:v1"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}

task MHC_parseBed {

    input {
        File bed
        Int spos = 0
        String mhc_path_name
        String output_prefix

        Int? preemptible_tries
    }


    command <<<
        set -eo pipefail

        python - --bed_file ~{bed} \
                 --spos ~{spos} \
                 --chromosome ~{mhc_path_name} \
                 --output_file ~{output_prefix} \
                 <<-'EOF'
        import gzip
        import argparse


        def import_bed(bed_file, spos, chromosome):
            locus_list = []
            gene_list = []
            hg38_locus = []
            if bed_file.endswith("gz"):
                with gzip.open(bed_file, "r") as f:
                    for line in f:
                        itemlist = line.strip().split("\t")
                        locus = "%s:%d-%d" % (chromosome, int(itemlist[1])-spos, int(itemlist[2]) - spos)
                        hg38_l = "%s:%d-%d" % (itemlist[0], int(itemlist[1]), int(itemlist[2]))
                        locus_list.append(locus)
                        gene_list.append(itemlist[3])
                        hg38_locus.append(hg38_l)
            else:
                with open(bed_file, "r") as f:
                    for line in f:
                        itemlist = line.strip().split("\t")
                        locus = "%s:%d-%d" % (chromosome, int(itemlist[1]) - spos, int(itemlist[2]) - spos)
                        hg38_l = "%s:%d-%d" % (itemlist[0], int(itemlist[1]), int(itemlist[2]))
                        hg38_locus.append(hg38_l)
                        locus_list.append(locus)
                        gene_list.append(itemlist[3])
            return (locus_list, gene_list, hg38_locus)


        def write_file(content, output_file):
            with open(output_file, "w") as f:
                for item in content:
                    f.write(f"{item}\n")

        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--bed_file',
                                type=str,
                                required=True)

            parser.add_argument('--spos',
                                type=int,
                                required=True)

            parser.add_argument('--chromosome',
                                type=str,
                                required=True)

            parser.add_argument('--output_file',
                                type=str,
                                required=True)

            args = parser.parse_args()

            locus_list, gene_list, hg38_locus = import_bed(args.bed_file, args.spos, args.chromosome)
            write_file(locus_list, args.output_file + "_locus.txt")
            write_file(gene_list, args.output_file + "_gene.txt")
            write_file(hg38_locus, args.output_file + "_hg38_locus.txt")


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
        Array[String] reflocuslist= read_lines("~{output_prefix}_hg38_locus.txt")
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
    }

    input {
        File bam
        File bai
        String locus
        String prefix = "subset"
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

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/vg:1.69.0"
        cpu: 1
        memory: "4GB"
        disks: "local-disk " + 100 + " HDD"
        preemptible: 0
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

task liftover_vcf {
    input {
        File vcf
        String current_contig
        String target_contig
        Int path_spos
        String prefix
    }

    command <<<
        set -euxo pipefail
        echo -e "~{current_contig}\t~{target_contig}" > chrom.map
        bcftools annotate --rename-chrs chrom.map ~{vcf} | awk -v spos=~{path_spos} 'BEGIN{OFS="\t"} /^#/ {print $0; next} {print $1, $2 + spos, $3, $4, $5, $6, $7, $8, $9, $10}' > test.vcf
        bcftools view test.vcf -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz


    >>>

    output {
        File lifted_vcf = "~{prefix}.vcf.gz"
        File lifted_tbi = "~{prefix}.vcf.gz.tbi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/hifiasm:0.25.0"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }
}
