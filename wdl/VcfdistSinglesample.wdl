version 1.0

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

struct OverlapMetricsOutputs {
    File inconsistent_tsv
    File metrics_tsv
}

workflow VcfdistAndSwitchEvaluation {
    input {
        String sample
        File truth_vcf
        File confident_regions_bed_file
        File eval_vcf
        File? eval_vcf_idx
        File gene_bed
        File reference_fasta
        File reference_fasta_fai

        File vcfdist_bed_file
        String? vcfdist_extra_args
        Int? vcfdist_mem_gb

        String? switch_filter_args
        String summarize_evaluations_docker
    }

    call parseBed {
        input:
            bed = gene_bed,
            output_prefix = sample
    }

    scatter (pair in zip(parseBed.locuslist, parseBed.genelist)) {
        String locus = pair.left
        String gene_name = pair.right

        call SubsetSampleFromVcf as SubsetSampleFromVcfEval { input:
            vcf = eval_vcf,
            vcf_idx = eval_vcf_idx,
            original_sample_name = "syndip",
            sample = sample,
            region = locus,
            bed_file = confident_regions_bed_file,
            reference_fasta_fai = reference_fasta_fai
        }

        call SubsetSampleFromVcf as SubsetSampleFromVcfTruth { input:
            vcf = truth_vcf,
            original_sample_name = "syndip",
            sample = sample,
            region = locus,
            bed_file = confident_regions_bed_file,
            reference_fasta_fai = reference_fasta_fai
        }

        call switch { input: 
            truth_bcf = SubsetSampleFromVcfTruth.single_sample_vcf,
            truth_bcf_index = SubsetSampleFromVcfTruth.single_sample_vcf_tbi,
            test_bcf = SubsetSampleFromVcfEval.single_sample_vcf,
            test_bcf_index = SubsetSampleFromVcfEval.single_sample_vcf_tbi,
            region = locus,
            outputprefix = sample + "_" + gene_name,
            filter_args = switch_filter_args,
            num_threads = 4
        }
    
        call Vcfdist { input:
            sample = sample + "_" + gene_name,
            eval_vcf = SubsetSampleFromVcfEval.single_sample_vcf,
            truth_vcf = SubsetSampleFromVcfTruth.single_sample_vcf,
            bed_file = vcfdist_bed_file,
            reference_fasta = reference_fasta,
            extra_args = vcfdist_extra_args,
            mem_gb = vcfdist_mem_gb
        }


        call SummarizeEvaluations { input:
            labels_per_vcf = [gene_name],
            vcfdist_outputs_per_vcf_and_sample = [[Vcfdist.outputs]],
            docker = summarize_evaluations_docker
        }
    }

    output {
        # stratification x sample
        Array[VcfdistOutputs] vcfdist_summary = Vcfdist.outputs
        Array[File] evaluation_summary_tsv = SummarizeEvaluations.evaluation_summary_tsv
        # per sample
        Array[Array[File]] switch_output_files = switch.output_files
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


task SubsetSampleFromVcf {
    input {
        File vcf
        File? vcf_idx
        String original_sample_name
        String sample
        String region
        File? bed_file
        File reference_fasta_fai
    }

    command <<<
        set -euxo pipefail

        if ! ~{defined(vcf_idx)}; then
            bcftools index ~{vcf}
        fi

        # must use -T bed_file to intersect with -r region properly
        bcftools view ~{vcf} \
            -s ~{original_sample_name} \
            -r ~{region} \
            ~{"-T " + bed_file} \
            -Oz -o ~{sample}.subset.g.vcf.gz
        echo ~{sample} > sample_name.txt
        bcftools reheader ~{sample}.subset.g.vcf.gz \
            -s sample_name.txt \
            --fai ~{reference_fasta_fai} \
            -o ~{sample}.subset.reheadered.g.vcf.gz
        bcftools index -t ~{sample}.subset.reheadered.g.vcf.gz
    >>>
    
    output {
        File single_sample_vcf = "~{sample}.subset.reheadered.g.vcf.gz"
        File single_sample_vcf_tbi = "~{sample}.subset.reheadered.g.vcf.gz.tbi"
    }

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task Vcfdist {
    input {
        String sample
        File eval_vcf
        File truth_vcf
        File bed_file
        File reference_fasta
        String? extra_args
        Int verbosity = 1

        Int disk_size_gb = ceil(size(truth_vcf, "GiB") + 10)
        Int mem_gb = 32
        Int cpu = 4
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        vcfdist \
            ~{eval_vcf} \
            ~{truth_vcf} \
            ~{reference_fasta} \
            -b ~{bed_file} \
            -v ~{verbosity} \
            ~{extra_args}

        for tsv in $(ls *.tsv); do mv $tsv ~{sample}.$tsv; done
        mv summary.vcf ~{sample}.summary.vcf
    >>>

    output {
        VcfdistOutputs outputs = {
            "summary_vcf": "~{sample}.summary.vcf",
            "precision_recall_summary_tsv": "~{sample}.precision-recall-summary.tsv",
            "precision_recall_tsv": "~{sample}.precision-recall.tsv",
            "query_tsv": "~{sample}.query.tsv",
            "truth_tsv": "~{sample}.truth.tsv",
            "phasing_summary_tsv": "~{sample}.phasing-summary.tsv",
            "switchflips_tsv": "~{sample}.switchflips.tsv",
            "superclusters_tsv": "~{sample}.superclusters.tsv",
            "phase_blocks_tsv": "~{sample}.phase-blocks.tsv"
        }
    }

    runtime {
        docker: "timd1/vcfdist:v2.5.3"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}

# see description at https://github.com/broadinstitute/lrma-aou1-panel-creation/pull/22#issue-2456295093
task CalculateOverlapMetrics {
    input {
        File vcf
        File vcf_idx
        String region
        String? phase_tag    # e.g., "PG" for kanpig integrated, "PS" for HiPhase, "NONE" for Shapeit4

        # docker needs bcftools, scikit-allel, and pandas
        String docker
        Int disk_size_gb = ceil(size(vcf, "GiB") + 10)
        Int mem_gb = 32
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        # split to biallelic
        bcftools norm \
            -r ~{region} \
            -m- --do-not-normalize ~{vcf} \
            -Oz -o split.vcf.gz
        bcftools index -t split.vcf.gz

        # join and subset to multiallelic
        bcftools norm \
            -r ~{region} \
            -m+ --do-not-normalize ~{vcf} | \
            bcftools view --min-alleles 3 \
                -Oz -o multi.vcf.gz
        bcftools index -t multi.vcf.gz

        python - --split_vcf split.vcf.gz \
                 --joined_multi_vcf multi.vcf.gz \
                 --phase_tag ~{phase_tag} \
                 <<-'EOF'
        import argparse
        import allel
        import numpy as np
        import pandas as pd
        from collections import Counter

        def calculate_metrics(split_vcf, joined_multi_vcf, phase_tag):
            joined_multi_callset = allel.read_vcf(joined_multi_vcf)
            split_callset = allel.read_vcf(split_vcf,
                                           fields=['samples', 'calldata/GT', 'CHROM', 'POS'] + ([] if phase_tag == 'NONE' else [f'calldata/{phase_tag}']))

            # V indexes all split records, v indexes only those that are in multiallelic sites
            is_multi_V = np.isin(split_callset['variants/POS'], joined_multi_callset['variants/POS'])
            gt_vsp = split_callset['calldata/GT'][is_multi_V]
            is_phased_vs = np.ones(gt_vsp.shape[:2]) if phase_tag == 'NONE' else split_callset[f'calldata/{phase_tag}'][is_multi_V] > 0
            chrom_v = split_callset['variants/CHROM'][is_multi_V]
            pos_v = split_callset['variants/POS'][is_multi_V]
            samples_s = split_callset['samples']

            # group split records by multiallelic site (indexed by m)
            m_splits = np.unique(pos_v, return_index=True)[1][1:]
            gt_mvsp = np.split(gt_vsp, m_splits)
            is_phased_mvs = np.split(is_phased_vs, m_splits)
            chrom_mv = np.split(chrom_v, m_splits)
            pos_mv = np.split(pos_v, m_splits)
            
            # calculate metrics and record inconsistent positions/samples
            metrics = Counter({'NUM_CONSISTENT_ALLELES': 0,
                               'NUM_INCONSISTENT_ALLELES': 0,
                               'NUM_CONSISTENT_SITES': 0,
                               'NUM_INCONSISTENT_SITES': 0})
            inconsistent_records = []

            for m, site_gt_vsp in enumerate(gt_mvsp):
                # site is considered inconsistent if at least one sample has phased genotype calls
                # that yield multiple alternate alleles on the same haplotype
                is_phased_alt_vsp = ((site_gt_vsp != 0) & (site_gt_vsp != -1)) * is_phased_mvs[m][:, :, np.newaxis]
                is_inconsistent_sp = (is_phased_alt_vsp).sum(axis=0) > 1
                is_inconsistent = is_inconsistent_sp.any()
                
                if is_inconsistent:
                    # find inconsistent alleles (i.e., split records indexed by v)
                    is_inconsistent_v = is_phased_alt_vsp[:, is_inconsistent_sp].any(axis=1)

                    # record genomic positions, samples, allele indices
                    chrom = chrom_mv[m][0]
                    start = pos_mv[m][0]
                    end = pos_mv[m][-1]
                    is_inconsistent_s = is_inconsistent_sp.any(axis=1)
                    samples = ','.join(samples_s[is_inconsistent_s])
                    allele_indices = ','.join(np.where(is_inconsistent_v)[0].astype(str))
                    inconsistent_records.append([chrom, start, end, samples, allele_indices])

                    # update counts
                    num_inconsistent_alleles = is_inconsistent_v.sum()
                    metrics['NUM_INCONSISTENT_ALLELES'] += num_inconsistent_alleles
                    metrics['NUM_CONSISTENT_ALLELES'] += len(site_gt_vsp) - num_inconsistent_alleles
                    metrics['NUM_INCONSISTENT_SITES'] += 1
                else:
                    metrics['NUM_CONSISTENT_ALLELES'] += len(site_gt_vsp)
                    metrics['NUM_CONSISTENT_SITES'] += 1

            # write inconsistent records
            inconsistent_df = pd.DataFrame(inconsistent_records,
                                           columns=['CHROM', 'START', 'END', 'SAMPLES', 'ALLELE_INDICES'])
            inconsistent_df.to_csv('inconsistent.tsv', sep='\t', index=False)

            # write in/consistent counts
            metrics_df = pd.DataFrame.from_records([dict(metrics)])
            metrics_df.to_csv('metrics.tsv', sep='\t', index=False,
                              columns=['NUM_INCONSISTENT_ALLELES', 'NUM_CONSISTENT_ALLELES',
                                       'NUM_INCONSISTENT_SITES', 'NUM_CONSISTENT_SITES'])


        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--split_vcf',
                                type=str)

            parser.add_argument('--joined_multi_vcf',
                                type=str)

            parser.add_argument('--phase_tag',
                                type=str)

            args = parser.parse_args()

            calculate_metrics(args.split_vcf,
                              args.joined_multi_vcf,
                              args.phase_tag)

        if __name__ == '__main__':
            main()
        EOF
    >>>

    output {
        OverlapMetricsOutputs outputs = {
            "inconsistent_tsv": "inconsistent.tsv",
            "metrics_tsv": "metrics.tsv"
        }
    }

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}

task SummarizeEvaluations {
    input {
        Array[String] labels_per_vcf
        Array[Array[VcfdistOutputs]] vcfdist_outputs_per_vcf_and_sample

        String docker
        Int disk_size_gb = 50
        Int boot_disk_size_gb = 15
        Int mem_gb = 8
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        python - --labels_per_vcf_txt ~{write_lines(labels_per_vcf)} \
                 --vcfdist_outputs_per_vcf_and_sample_json ~{write_json(vcfdist_outputs_per_vcf_and_sample)} \
                 <<-'EOF'
        import argparse
        import json
        import pandas as pd

        def summarize(labels_per_vcf_txt,
                      vcfdist_outputs_per_vcf_and_sample_json,
                      overlap_metrics_outputs_per_vcf):
            with open(labels_per_vcf_txt) as f:
                labels = f.read().splitlines()

            with open(vcfdist_outputs_per_vcf_and_sample_json) as f:
                vcfdist_outputs_per_vcf_and_sample = json.load(f)

            summary_dict = {}
            for i, label in enumerate(labels):
                summary_dict[label] = {}
                summary_dict[label]['NUM_VCFDIST_SAMPLES'] = len(vcfdist_outputs_per_vcf_and_sample[i])
                summary_dict[label].update(summarize_vcfdist_outputs_over_samples(vcfdist_outputs_per_vcf_and_sample[i]))
                if overlap_metrics_outputs_per_vcf:
                    summary_dict[label].update(summarize_overlap_metrics_outputs(overlap_metrics_outputs_per_vcf[i]))

            pd.DataFrame.from_dict(summary_dict, orient='index').to_csv('evaluation_summary.tsv', sep='\t', float_format="%.4f")

        def summarize_vcfdist_outputs_over_samples(vcfdist_outputs_per_sample):
            precision_recall_metrics_dict = {}
            for s, vcfdist_outputs in enumerate(vcfdist_outputs_per_sample):
                precision_recall_metrics_dict[s] = {}
                pr_metrics_df = pd.read_csv(vcfdist_outputs['precision_recall_summary_tsv'], sep='\t', index_col=[0, 1])
                for var_type in ['SNP', 'INDEL', 'SV']:
                    var_type_metrics_dict = pr_metrics_df.loc[var_type, 'NONE'][['TRUTH_TP', 'QUERY_TP', 'TRUTH_FN', 'QUERY_FP', 'PREC', 'RECALL', 'F1_SCORE']].add_prefix(f'{var_type}_').to_dict()
                    precision_recall_metrics_dict[s].update(var_type_metrics_dict)
            return pd.DataFrame.from_dict(precision_recall_metrics_dict, orient='index').mean(axis=0)


        def summarize_overlap_metrics_outputs(overlap_metrics_outputs):
            return pd.read_csv(overlap_metrics_outputs['metrics_tsv'], sep='\t').iloc[0].to_dict()

        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--labels_per_vcf_txt',
                                type=str)

            parser.add_argument('--vcfdist_outputs_per_vcf_and_sample_json',
                                type=str)

            parser.add_argument('--overlap_metrics_outputs_per_vcf_json',
                                type=str)

            args = parser.parse_args()

            summarize(args.labels_per_vcf_txt,
                      args.vcfdist_outputs_per_vcf_and_sample_json,
                      args.overlap_metrics_outputs_per_vcf_json)

        if __name__ == '__main__':
            main()
        EOF
    >>>

    output {
        File evaluation_summary_tsv = "evaluation_summary.tsv"
    }

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        bootDiskSizeGb: boot_disk_size_gb
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}

task switch {
    input{
        File truth_bcf
        File truth_bcf_index
        File test_bcf
        File test_bcf_index
        String region
        String outputprefix
        Int num_threads
        String? filter_args = "-i 'abs(strlen(ALT)-strlen(REF))<5'"
    }
    command <<<
        set -euxo pipefail
        bcftools view ~{filter_args} ~{truth_bcf} -Oz -o truth.short_filtered.bcf
        bcftools index truth.short_filtered.bcf
        bcftools view ~{filter_args} ~{test_bcf} -Oz -o test.short_filtered.bcf
        bcftools index test.short_filtered.bcf

        switch_static --validation truth.short_filtered.bcf --estimation test.short_filtered.bcf --region ~{region} --output ~{outputprefix} --thread ~{num_threads}
    >>>

    output{
        Array[File] output_files = glob("*")
    }

    Int disk_size = 100 + ceil(2 * (size(truth_bcf, "GiB") + size(test_bcf, "GiB")))

    runtime {
        cpu: 4
        memory: "16 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/shapeit5:v1"
    }
}