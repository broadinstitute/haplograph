
version 1.0

workflow VG_remap_reads_hifiasm {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File gfa_file
        Array[Int] desiredCoverages
        Int hifiasm_mem
        Int hifiasm_thread
        File reference_fa
        String prefix
        Int giraffe_threads
        String sample_id
        String giraffe_extra_arg
        String graph_ref_path
        String gene_name = "mhc"
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
            whole_genome_bam = whole_genome_bam,
            whole_genome_bai = whole_genome_bai,
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

    call CalculateCoverage {
        input:
            bam = vg_graph_surject.bam_file,
            bai = vg_graph_surject.bam_index_file,
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

        call hifiasm_asm {
            input:
                bam = downsampleBam.downsampled_bam,
                bai = downsampleBam.downsampled_bai,
                prefix = prefix + "_" + gene_name + "_" + desiredCoverage,
                num_cpus = hifiasm_thread,
                mem_gb = hifiasm_mem
        }


    }
  
    output {
        File graph_bam = vg_graph_surject.bam_file
        File graph_bai = vg_graph_surject.bam_index_file
        Array[File] hifiasm_asm_hap1 = hifiasm_asm.assembly_hap1
        Array[File] hifiasm_asm_hap2 = hifiasm_asm.assembly_hap2
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
        File whole_genome_bam
        File whole_genome_bai
        Int thread_num
        String read_type = "hifi"
        String in_sample_name
        String extra_args
        String out_prefix
    }
    Int half_cores = thread_num / 2

    command <<<
        set -euxo pipefail

        samtools collate -Oun128 ~{whole_genome_bam} | samtools fastq -OT RG,BC,Mm,Ml - \
        | /vg/vg giraffe --parameter-preset ~{read_type} \
                       --progress --track-provenance \
                       -t ~{thread_num} \
                       -Z ~{graph_index} \
                       -d ~{graph_dist} \
                       -m ~{graph_min} \
                       -z ~{graph_zipcode} \
                       -f - \
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
        disks: "local-disk " + 1000 + " HDD"
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
        prefix: "prefix for output bam and bai file names"
    }

    input {
        File bam
        File bai
        String prefix = "subset"
    }



    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} > ~{prefix}.bam
        samtools index ~{prefix}.bam
        samtools depth ~{prefix}.bam | awk '{sum+=$3} END {print sum/NR}' > coverage.txt

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
