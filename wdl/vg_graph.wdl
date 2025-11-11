
version 1.0

workflow VG_extract_bam {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File gfa_file
        String prefix
        String locus
        String vg_locus
        Int giraffe_threads
        String sample_id
        String giraffe_extra_arg
        String graph_ref_path
    }

    call SubsetBam {
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
            locus = locus,
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

    call haplograph {
        input:
            bam = vg_graph_surject.bam_file,
            bai = vg_graph_surject.bam_index_file,
            reference_fa = vg_graph_index.graph_reference_fa,
            prefix = prefix,
            locus = vg_locus
    }
    
    
    output {
        File graph_bam = vg_graph_surject.bam_file
        File graph_bai = vg_graph_surject.bam_index_file
        File vcf = haplograph.germline_vcf_file

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
        String? extra_args
        Int verbosity = 1

        Int disk_size_gb = ceil(size(truth_vcf, "GiB") + 10)
        Int mem_gb = 32
        Int cpu = 8
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
