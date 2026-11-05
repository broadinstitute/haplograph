version 1.0

workflow BuildSparseGRM {
    meta {
        description: "Filter per-chromosome VCFs to a sparse biallelic SNP panel, subset to an analysis sample list, and build a PLINK bed dataset plus a SAIGE sparse GRM. Long alleles are removed before plink2 --set-all-var-ids, which aborts when any allele code exceeds 23 characters. Pass chrom_vcfs in concatenation order (chr1, chr2, ...)."
    }

    input {
        # Sample IDs to keep. Matched against the VCF sample column (PLINK IID after --double-id).
        File samplelist
        # One VCF per chromosome. Filtered, then concatenated in this order before PLINK import.
        Array[File] chrom_vcfs
        # Optional indexes, parallel to chrom_vcfs (.tbi or .csi). Built in the task when omitted.
        Array[File]? chrom_vcf_indices

        Float relatedness_cutoff = 0.05
        Int num_random_markers = 2000

        # Applied per chromosome with bcftools before plink2 import.
        # plink2 checks allele length while writing the temporary .psam,
        # before --maf / --snps-only, so the VCF itself has to be filtered.
        # Biallelic SNPs stay under the 23-character --set-all-var-ids limit.
        Float min_maf = 0.01
        Float max_missing = 0.15
        # LD prune on the analysis samples after import. Set ld_prune=false to
        # keep every SNP that passed the VCF filter.
        Boolean ld_prune = true
        Int ld_window = 50
        Int ld_step = 5
        Float ld_r2 = 0.2

        # Image must provide bcftools, plink2, Rscript, and createSparseGRM.R
        # (/usr/local/bin or /opt/SAIGE/extdata).
        String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
        Int cpu = 8
        Int memory_gb = 32
        # Each chromosome is filtered and the raw VCF deleted before the next
        # one is read, so peak is the localized inputs plus one filtered copy
        # rather than a full unfiltered concat.
        Int disk_gb = ceil(size(chrom_vcfs, "GB") * 2) + 100
        Int preemptible = 0
    }

    call BuildSparseGRMAndPlink {
        input:
            analysis_samples = samplelist,
            chrom_vcfs = chrom_vcfs,
            chrom_vcf_indices = chrom_vcf_indices,
            relatedness_cutoff = relatedness_cutoff,
            num_random_markers = num_random_markers,
            min_maf = min_maf,
            max_missing = max_missing,
            ld_prune = ld_prune,
            ld_window = ld_window,
            ld_step = ld_step,
            ld_r2 = ld_r2,
            docker = docker,
            cpu = cpu,
            memory_gb = memory_gb,
            disk_gb = disk_gb,
            preemptible = preemptible
    }

    output {
        File plink_bed = BuildSparseGRMAndPlink.plink_bed
        File plink_bim = BuildSparseGRMAndPlink.plink_bim
        File plink_fam = BuildSparseGRMAndPlink.plink_fam
        File sparseGRM = BuildSparseGRMAndPlink.sparseGRM
        File sparseGRM_IDlist = BuildSparseGRMAndPlink.sparseGRM_IDlist
        File samples_kept = BuildSparseGRMAndPlink.samples_kept
        File paths_tsv = BuildSparseGRMAndPlink.paths_tsv
    }
}

task BuildSparseGRMAndPlink {
    input {
        File analysis_samples
        Array[File] chrom_vcfs
        Array[File]? chrom_vcf_indices
        Float relatedness_cutoff = 0.05
        Int num_random_markers = 2000
        Float min_maf = 0.01
        Float max_missing = 0.15
        Boolean ld_prune = true
        Int ld_window = 50
        Int ld_step = 5
        Float ld_r2 = 0.2
        String docker
        Int cpu = 8
        Int memory_gb = 32
        Int disk_gb = 400
        Int preemptible = 0
    }

    Array[File] vcf_indices = select_first([chrom_vcf_indices, []])

    command <<<
        set -euo pipefail

        # R, bcftools, and plink2 spill into /tmp by default. On Terra that
        # disk is the small boot disk.
        export TMPDIR="${PWD}/tmp" TMP="${PWD}/tmp" TEMP="${PWD}/tmp"
        mkdir -p "${TMPDIR}" grm_build grm_build/filt

        say_disk() {
            echo "[disk] $1: $(df -Pk . | awk 'NR==2{printf "%.1f GiB used, %.1f GiB free", $3/1048576, $4/1048576}')"
        }

        if [[ -x /usr/local/bin/createSparseGRM.R ]]; then
            CREATE_SPARSE_GRM_R=/usr/local/bin/createSparseGRM.R
        elif [[ -f /opt/SAIGE/extdata/createSparseGRM.R ]]; then
            CREATE_SPARSE_GRM_R=/opt/SAIGE/extdata/createSparseGRM.R
        elif command -v createSparseGRM.R >/dev/null 2>&1; then
            CREATE_SPARSE_GRM_R=$(command -v createSparseGRM.R)
        else
            echo "createSparseGRM.R not found in /usr/local/bin or /opt/SAIGE/extdata" >&2
            exit 1
        fi

        N_THREADS=~{cpu}
        say_disk "start"

        # One ID per line. Drop blanks and duplicate IDs, keep first-seen order.
        awk 'NF && !seen[$1]++ { print $1 }' "~{analysis_samples}" > grm_build/analysis_samples.txt
        N_WANT=$(wc -l < grm_build/analysis_samples.txt | tr -d ' ')
        if [[ "${N_WANT}" -eq 0 ]]; then
            echo "samples is empty" >&2
            exit 1
        fi

        VCF_LIST="~{write_lines(chrom_vcfs)}"
        IDX_LIST="~{write_lines(vcf_indices)}"
        N_VCF=$(grep -c . "${VCF_LIST}" || true)
        if [[ "${N_VCF}" -lt 1 ]]; then
            echo "chrom_vcfs is empty" >&2
            exit 1
        fi
        N_IDX=0
        if [[ -s "${IDX_LIST}" ]]; then
            N_IDX=$(grep -c . "${IDX_LIST}" || true)
            if [[ "${N_IDX}" -ne "${N_VCF}" ]]; then
                echo "chrom_vcf_indices (${N_IDX}) must be parallel to chrom_vcfs (${N_VCF})" >&2
                exit 1
            fi
        fi

        # plink2 applies --set-all-var-ids while writing
        # grm_build/all-temporary.psam, before --maf or --snps-only. Any allele
        # longer than 23 characters aborts that step ("allele codes too long").
        # Filter the VCF first: analysis samples, biallelic SNPs, MAF, and
        # missingness. One chromosome at a time, deleting the raw VCF as soon
        # as its filtered copy exists.
        MIN_MAF="~{min_maf}"
        MAX_MISSING="~{max_missing}"
        FILTER_EXPR="MAF>=${MIN_MAF} && F_MISSING<=${MAX_MISSING} && strlen(REF)==1 && strlen(ALT)==1"
        echo "VCF filter: biallelic SNPs, ${FILTER_EXPR}"

        bcftools query -l "$(head -n 1 "${VCF_LIST}")" > grm_build/vcf_samples.txt
        awk 'NR==FNR { want[$1]=1; next } ($1 in want) { print $1 }' \
            grm_build/analysis_samples.txt grm_build/vcf_samples.txt \
            > grm_build/samples_in_vcf.txt
        N_IN_VCF=$(wc -l < grm_build/samples_in_vcf.txt | tr -d ' ')
        echo "Analysis samples present in the VCF header: ${N_IN_VCF} / ${N_WANT}"
        if [[ "${N_IN_VCF}" -eq 0 ]]; then
            echo "No analysis samples are present in the VCFs" >&2
            exit 1
        fi
        if [[ "${N_IN_VCF}" -lt "${N_WANT}" ]]; then
            echo "WARNING: $((N_WANT - N_IN_VCF)) analysis samples are absent from the VCF header and will be dropped" >&2
            awk 'NR==FNR { have[$1]=1; next } !($1 in have) { print $1 }' \
                grm_build/vcf_samples.txt grm_build/analysis_samples.txt \
                | head -n 20 >&2 || true
        fi

        FILTERED_LIST=grm_build/filtered_vcfs.txt
        : > "${FILTERED_LIST}"
        N_VAR_TOTAL=0
        i=0
        while IFS= read -r vcf || [[ -n "${vcf}" ]]; do
            [[ -z "${vcf}" ]] && continue
            i=$((i + 1))
            filt="grm_build/filt/${i}.vcf.gz"
            echo "Filtering chromosome file ${i}/${N_VCF}: ${vcf}"
            bcftools view \
                -S grm_build/samples_in_vcf.txt \
                -v snps -m2 -M2 \
                -i "${FILTER_EXPR}" \
                -Oz -o "${filt}" \
                --threads "${N_THREADS}" \
                "${vcf}"
            bcftools index -t "${filt}"
            n_kept=$(bcftools index -n "${filt}" | tr -d ' ')
            echo "  kept ${n_kept} variants"
            N_VAR_TOTAL=$((N_VAR_TOTAL + n_kept))
            echo "${filt}" >> "${FILTERED_LIST}"
            rm -f "${vcf}" "${vcf}.tbi" "${vcf}.csi" || true
            if [[ "${N_IDX}" -gt 0 ]]; then
                idx=$(sed -n "${i}p" "${IDX_LIST}")
                [[ -n "${idx}" ]] && rm -f "${idx}" || true
            fi
            say_disk "after filtering chromosome ${i}"
        done < "${VCF_LIST}"

        echo "Variants after VCF filter: ${N_VAR_TOTAL}"
        if [[ "${N_VAR_TOTAL}" -lt "~{num_random_markers}" ]]; then
            echo "Only ${N_VAR_TOTAL} variants survived filtering; createSparseGRM needs at least ~{num_random_markers}" >&2
            exit 1
        fi

        if [[ "${N_VCF}" -eq 1 ]]; then
            INPUT_VCF=$(head -n 1 "${FILTERED_LIST}")
        else
            bcftools concat -f "${FILTERED_LIST}" -Oz -o grm_build/concat.vcf.gz --threads "${N_THREADS}"
            bcftools index -t grm_build/concat.vcf.gz
            INPUT_VCF=grm_build/concat.vcf.gz
            rm -rf grm_build/filt
            say_disk "after concat of filtered VCFs"
        fi

        # FLARE ancestry VCFs often declare ##FORMAT=<ID=GT,...> twice, and
        # header merge during concat can keep both. PLINK2 rejects that.
        # Skip the reheader copy when there is a single GT line.
        N_GT=$(bcftools view -h "${INPUT_VCF}" | grep -c '^##FORMAT=<ID=GT,' || true)
        echo "FORMAT/GT header lines: ${N_GT}"
        if [[ "${N_GT}" -gt 1 ]]; then
            bcftools view -h "${INPUT_VCF}" \
                | awk '/^##FORMAT=<ID=GT,/{ if (gt++) next } { print }' \
                > grm_build/header.fixed.txt
            bcftools reheader -h grm_build/header.fixed.txt -o grm_build/plink_in.vcf.gz "${INPUT_VCF}"
            bcftools index -t grm_build/plink_in.vcf.gz \
                || bcftools index -c grm_build/plink_in.vcf.gz
            rm -f grm_build/concat.vcf.gz grm_build/concat.vcf.gz.tbi grm_build/concat.vcf.gz.csi
            rm -rf grm_build/filt
            INPUT_VCF=grm_build/plink_in.vcf.gz
            say_disk "after reheader"
        else
            echo "Single GT header line: skipping the reheader copy."
        fi

        plink2 --vcf "${INPUT_VCF}" \
            --double-id \
            --make-bed \
            --set-all-var-ids '@:#:$r:$a' \
            --out grm_build/all \
            --threads "${N_THREADS}"
        say_disk "after plink2 import"
        rm -rf grm_build/filt
        rm -f grm_build/concat.vcf.gz grm_build/concat.vcf.gz.tbi grm_build/concat.vcf.gz.csi
        rm -f grm_build/plink_in.vcf.gz grm_build/plink_in.vcf.gz.tbi grm_build/plink_in.vcf.gz.csi
        say_disk "after freeing VCFs"

        # Keep the intersection of the requested IDs and the imported FAM.
        # --double-id sets FID and IID to the VCF sample name.
        awk 'NR==FNR { have[$2]=1; next } ($1 in have) { print $1 }' \
            grm_build/all.fam grm_build/analysis_samples.txt \
            > grm_build/analysis_samples.intersect.txt
        awk 'NR==FNR { want[$1]=1; next } ($2 in want) { print $1, $2 }' \
            grm_build/analysis_samples.txt grm_build/all.fam \
            > grm_build/keep_iid.txt
        N_KEEP=$(wc -l < grm_build/keep_iid.txt | tr -d ' ')
        echo "Keeping ${N_KEEP} / ${N_WANT} analysis samples"
        if [[ "${N_KEEP}" -eq 0 ]]; then
            echo "No analysis samples are present in the VCFs" >&2
            exit 1
        fi
        if [[ "${N_KEEP}" -lt "${N_WANT}" ]]; then
            echo "WARNING: $((N_WANT - N_KEEP)) analysis samples are absent from the VCFs and were dropped" >&2
            awk 'NR==FNR { have[$2]=1; next } !($1 in have) { print $1 }' \
                grm_build/all.fam grm_build/analysis_samples.txt \
                | head -n 20 >&2 || true
        fi

        mkdir -p saige_grm
        plink2 --bfile grm_build/all \
            --keep grm_build/keep_iid.txt \
            --make-bed \
            --out saige_grm/plink \
            --threads "${N_THREADS}"
        rm -f grm_build/all.bed grm_build/all.bim grm_build/all.fam
        N_MARKERS=$(wc -l < saige_grm/plink.bim | tr -d ' ')
        if [[ "~{ld_prune}" == "true" ]]; then
            echo "LD pruning ${N_MARKERS} markers: window=~{ld_window} step=~{ld_step} r2=~{ld_r2}"
            plink2 --bfile saige_grm/plink \
                --indep-pairwise ~{ld_window} ~{ld_step} ~{ld_r2} \
                --out grm_build/ldprune \
                --threads "${N_THREADS}"
            N_PRUNED=$(wc -l < grm_build/ldprune.prune.in | tr -d ' ')
            echo "LD prune kept ${N_PRUNED} / ${N_MARKERS} markers"
            if [[ "${N_PRUNED}" -lt "~{num_random_markers}" ]]; then
                echo "LD prune left ${N_PRUNED} markers; createSparseGRM needs at least ~{num_random_markers}" >&2
                exit 1
            fi
            plink2 --bfile saige_grm/plink \
                --extract grm_build/ldprune.prune.in \
                --make-bed \
                --out saige_grm/plink.pruned \
                --threads "${N_THREADS}"
            rm -f saige_grm/plink.bed saige_grm/plink.bim saige_grm/plink.fam
            mv saige_grm/plink.pruned.bed saige_grm/plink.bed
            mv saige_grm/plink.pruned.bim saige_grm/plink.bim
            mv saige_grm/plink.pruned.fam saige_grm/plink.fam
            N_MARKERS="${N_PRUNED}"
        fi
        say_disk "before createSparseGRM"
        echo "Markers entering createSparseGRM: ${N_MARKERS}"

        Rscript "${CREATE_SPARSE_GRM_R}" \
            --plinkFile=saige_grm/plink \
            --nThreads="${N_THREADS}" \
            --outputPrefix=saige_grm/sparseGRM \
            --numRandomMarkerforSparseKin=~{num_random_markers} \
            --relatednessCutoff=~{relatedness_cutoff} \
            --minMAFforGRM=~{min_maf} \
            --maxMissingRateforGRM=~{max_missing}

        shopt -s nullglob
        mtx_hits=(saige_grm/sparseGRM*.sparseGRM.mtx)
        id_hits=(saige_grm/sparseGRM*.sparseGRM.mtx.sampleIDs.txt)
        shopt -u nullglob
        SPARSE_MTX="${mtx_hits[0]:-}"
        SPARSE_IDS="${id_hits[0]:-}"
        if [[ -z "${SPARSE_MTX}" || -z "${SPARSE_IDS}" ]]; then
            echo "Failed to locate createSparseGRM outputs under saige_grm" >&2
            ls -la saige_grm >&2 || true
            exit 1
        fi

        # Materialize real files. SAIGE embeds the cutoff and marker count in
        # the filename; stable names are what downstream SAIGE tasks take.
        cp -L "${SPARSE_MTX}" sparseGRM.mtx
        cp -L "${SPARSE_IDS}" sparseGRM.sampleIDs.txt
        cp saige_grm/plink.bed plink.bed
        cp saige_grm/plink.bim plink.bim
        cp saige_grm/plink.fam plink.fam
        cp grm_build/analysis_samples.intersect.txt samples_kept.txt

        {
            printf 'key\tpath\n'
            printf 'plink_prefix\t%s\n' saige_grm/plink
            printf 'sparse_grm_mtx\t%s\n' "${SPARSE_MTX}"
            printf 'sparse_grm_sample_ids\t%s\n' "${SPARSE_IDS}"
            printf 'relatedness_cutoff\t%s\n' "~{relatedness_cutoff}"
            printf 'num_random_markers\t%s\n' "~{num_random_markers}"
            printf 'min_maf\t%s\n' "~{min_maf}"
            printf 'max_missing\t%s\n' "~{max_missing}"
            printf 'ld_prune\t%s\n' "~{ld_prune}"
            printf 'ld_window\t%s\n' "~{ld_window}"
            printf 'ld_step\t%s\n' "~{ld_step}"
            printf 'ld_r2\t%s\n' "~{ld_r2}"
            printf 'n_markers\t%s\n' "${N_MARKERS}"
            printf 'n_samples_requested\t%s\n' "${N_WANT}"
            printf 'n_samples_kept\t%s\n' "${N_KEEP}"
        } > saige_sparse_grm.paths.tsv

        say_disk "done"
        echo "Wrote PLINK plink.bed/.bim/.fam and sparse GRM:"
        echo "  ${SPARSE_MTX}"
        echo "  ${SPARSE_IDS}"
    >>>

    output {
        File plink_bed = "plink.bed"
        File plink_bim = "plink.bim"
        File plink_fam = "plink.fam"
        File sparseGRM = "sparseGRM.mtx"
        File sparseGRM_IDlist = "sparseGRM.sampleIDs.txt"
        File samples_kept = "samples_kept.txt"
        File paths_tsv = "saige_sparse_grm.paths.tsv"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}
