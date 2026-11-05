version 1.0

workflow SAIGE_phewas_null {
    meta {
        description: "PheWAS workflow for variant association using SAIGE — each task processes a batch of phecodes."
    }

    input {
        Array[String] phecode_list
        File homoplasmic_vcf
        File homoplasmic_vcf_csi
        File? heteroplasmic_vcf
        File? heteroplasmic_vcf_csi
        File? GroupFile
        File sparseGRM
        File sparseGRM_IDlist
        Array[File] null_model_rdas_list
        Array[File]? variance_ratio_txts_list
        Array[String] chromosome_list
        # Step 1 does not run here; the null models already encode the covariates
        # and the sample set (including the one-sex restriction a sex-specific
        # phecode needs). covariate_list / categorical_covariate_list /
        # phenotype_file / trait_type / step1_memory are accepted only so an
        # input JSON written for SAIGE_phewas.wdl also validates against this
        # workflow. Changing them here changes nothing — refit in SAIGE_phewas.wdl.
        Array[String] covariate_list  = ["PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","sex","age","age2","age_sex","age2_sex"]
        Array[String] categorical_covariate_list = ["sex"]
        File phenotype_file
        String trait_type
        String output_prefix
        String saige_docker
        String vcffield
        Boolean is_overridefilp
        # Firth effect sizes help rare binary traits but can crash SAIGE
        # (Mat::operator index out of bounds / glm separation). Default TRUE
        # matches prior behavior; set FALSE to skip Firth when Step 2 dies on
        # a subset of phecodes.
        Boolean is_Firth_beta = true

        # How Step 2 fills missing dosages (DS "."): best_guess = round(2*AF) (SAIGE's
        # default) turns a missing HLA/KIR call into 0 for any allele under 25% frequency;
        # mean = 2*AF, which drops the sample's contribution to the score test; minor = 0.
        # Has no effect on a VCF without missing values.
        String impute_method = "mean"

        String memory = "8G"
        # Step 1 (fitNULLGLMM) loads the full AoU sparse GRM + plink into R; 8G gets OOM-killed
        # (SIGKILL / exit 137), leaving 0-byte .rda outputs. Keep Step 1 on its own larger knob.
        String step1_memory = "32G"
        
        Float single_variant_minimal_af = 0.001
        Int single_variant_min_mac = 20
        Float single_variant_max_missing = 0.1
        Float gene_set_minimal_af = 0
        Float gene_set_min_mac = 0.5
        Float gene_set_max_missing = 0.5

        # MAC category bounds for categorical variance ratios, shared by Step 1 and Step 2 so the
        # two cannot drift. SAIGE requires length(min_exclude) == length(max_include) + 1, and
        # length(min_exclude) == the number of "null" rows Step 1 writes to the varianceRatio file.
        # If these disagree between steps, SAIGE applies each ratio to the wrong MAC bin silently.
        String cate_var_ratio_min_mac_exclude = "20,50"
        String cate_var_ratio_max_mac_include = "100"

        Int cpu = 2
        # 50G is tight: AoU sparse GRM + plink + 15 null .rda files. Empty Step 1
        # outputs are often a full /tmp or a full local-disk, not a SAIGE bug.
        String disk_size = "local-disk 200 HDD"
        Int preemptible = 1
    }
    
    # Step 2a: single-variant test (homoplasmic VCF) for all phecodes
    call RunStep2_singlevariant {
        input:
            vcf              = homoplasmic_vcf,
            vcf_csi          = homoplasmic_vcf_csi,
            variance_ratios  = variance_ratio_txts_list,
            GMMATmodelFiles  = null_model_rdas_list,
            sparseGRM        = sparseGRM,
            sparseGRM_IDlist = sparseGRM_IDlist,
            vcffield         = vcffield,
            chromo_list      = chromosome_list,
            phecode_list     = phecode_list,
            output_prefix    = output_prefix,
            minimal_af       = single_variant_minimal_af,
            min_mac          = single_variant_min_mac,
            max_missing      = single_variant_max_missing,
            is_overridefilp  = is_overridefilp,
            is_Firth_beta    = is_Firth_beta,
            impute_method    = impute_method,
            cate_var_ratio_min_mac_exclude = cate_var_ratio_min_mac_exclude,
            cate_var_ratio_max_mac_include = cate_var_ratio_max_mac_include,
            memory           = memory,
            saige_docker     = saige_docker,
            cpu              = cpu,
            disk             = disk_size,
            preemptible      = preemptible
    }

    # Step 2b: gene-set test (heteroplasmic VCF) for all phecodes
    if (defined(heteroplasmic_vcf)){
        call RunStep2_geneset {
            input:
                vcf              = select_first([heteroplasmic_vcf, ""]),
                vcf_csi          = select_first([heteroplasmic_vcf_csi, ""]),
                variance_ratios  = variance_ratio_txts_list,
                GMMATmodelFiles  = null_model_rdas_list,
                sparseGRM        = sparseGRM,
                sparseGRM_IDlist = sparseGRM_IDlist,
                vcffield         = vcffield,
                chromo_list      = chromosome_list,
                phecode_list     = phecode_list,
                output_prefix    = output_prefix,
                minimal_af       = gene_set_minimal_af,
                min_mac          = gene_set_min_mac,
                max_missing      = gene_set_max_missing,
                GroupFile        = select_first([GroupFile, ""]),
                is_Firth_beta    = is_Firth_beta,
                impute_method    = impute_method,
                cate_var_ratio_min_mac_exclude = cate_var_ratio_min_mac_exclude,
                cate_var_ratio_max_mac_include = cate_var_ratio_max_mac_include,
                memory           = memory,
                saige_docker     = saige_docker,
                cpu              = cpu,
                disk             = disk_size,
                preemptible      = preemptible
        }
    }


    output {
        Array[File] singlevariant_outputs  = RunStep2_singlevariant.outputfiles
        File        singlevariant_failed   = RunStep2_singlevariant.failed_list
        Array[File]? geneset_outputs        = RunStep2_geneset.outputfiles
        File?       geneset_failed          = RunStep2_geneset.failed_list
    }
}


# ---------------------------------------------------------------------------
# Task: step 2 single-variant SPA test for every phecode in phecode_list
# ---------------------------------------------------------------------------
task RunStep2_singlevariant {

    input {
        File vcf
        File vcf_csi
        Array[File]? variance_ratios
        Array[File] GMMATmodelFiles
        # Step 1 in this workflow always fits with --useSparseGRMtoFitNULL=TRUE, so Step 2 must be
        # given the same sparse GRM — otherwise SAIGE silently falls back to the "null" variance
        # ratios instead of the sparse-GRM ones (the guard for this is commented out upstream).
        File? sparseGRM
        File? sparseGRM_IDlist
        String vcffield
        Array[String] chromo_list
        Array[String] phecode_list
        String output_prefix
        Boolean is_overridefilp = false
        Boolean is_Firth_beta = true
        String impute_method = "mean"
        Float minimal_af
        Int min_mac
        Float max_missing
        String cate_var_ratio_min_mac_exclude = "10,20.5"
        String cate_var_ratio_max_mac_include = "20.5"
        # Runtime parameters
        String memory
        String saige_docker
        Int cpu
        String disk
        Int preemptible
    }
    # fastTest needs VR "null" (no-GRM) rows; without a VR file SAIGE cannot use it.
    Boolean is_fastTest = defined(variance_ratios)
    command <<<
        # No -e: one phecode/chrom SAIGE crash must not abort the batch (same
        # pattern as RunFitNullGLMM). Setup errors still exit explicitly below.
        set -uxo pipefail

        export TMPDIR="${PWD}/tmp" TMP="${PWD}/tmp" TEMP="${PWD}/tmp"
        mkdir -p "${TMPDIR}"
        : > failed_step2_singlevar.txt
        n_ok=0
        n_fail=0

        phecodes=(~{sep=' ' phecode_list})
        # Optional VR: undefined/empty → run without --varianceRatioFile (ratio=1).
        variance_ratios=(~{sep=' ' select_first([variance_ratios, []])})
        model_files=(~{sep=' ' GMMATmodelFiles})
        chromo_list=(~{sep=' ' chromo_list})
        have_vr=0
        if [[ ${#variance_ratios[@]} -gt 0 ]]; then
            have_vr=1
        else
            echo "WARNING: no variance-ratio files; Step 2 will omit --varianceRatioFile (is_fastTest=~{is_fastTest})" >&2
        fi

        # Pair phecode ↔ null model by file name, not by array index. The caller
        # normally passes null_model_rdas_list straight from a previous run's
        # glob() output, whose order is Cromwell's and is neither sorted nor the
        # order of a hand-written phecode_list. Zipping those by index labels
        # each result with another phecode's name, so a female-only null model
        # ends up scoring a both-sex phenotype (and vice versa) — which both
        # creates inflation and hides it.
        find_for() {
            local want="$1"; shift
            local f
            for f in "$@"; do
                if [[ "$(basename "${f}")" == *"${want}" && -s "${f}" ]]; then
                    printf '%s' "${f}"
                    return 0
                fi
            done
            return 1
        }

        if [[ ${#phecodes[@]} -eq 0 ]]; then
            echo "ERROR: phecode_list is empty" >&2
            exit 1
        fi
        if [[ ${#model_files[@]} -eq 0 ]]; then
            echo "ERROR: null_model_rdas_list is empty" >&2
            exit 1
        fi

        for i in "${!phecodes[@]}"; do
            phecode="${phecodes[$i]}"
            model_file="$(find_for "_step1Out_${phecode}.rda" "${model_files[@]}" || true)"
            variance_ratio=""
            if [[ -z "${model_file}" ]]; then
                # Pre-1.x outputs, or a prefix this task cannot reconstruct. Only
                # then fall back to the index, and say so loudly.
                if [[ ${#model_files[@]} -eq ${#phecodes[@]} ]]; then
                    model_file="${model_files[$i]}"
                    echo "WARNING: no .rda named *_step1Out_${phecode}.rda; falling back to position ${i} ($(basename "${model_file}")) — verify this pairing" >&2
                fi
            fi
            if [[ "${have_vr}" -eq 1 ]]; then
                variance_ratio="$(find_for "_step1Out_${phecode}.varianceRatio.txt" "${variance_ratios[@]}" || true)"
                if [[ -z "${variance_ratio}" && ${#variance_ratios[@]} -eq ${#phecodes[@]} ]]; then
                    variance_ratio="${variance_ratios[$i]}"
                    echo "WARNING: no varianceRatio named for ${phecode}; falling back to position ${i}" >&2
                fi
            fi
            if [[ -z "${model_file}" || ! -s "${model_file}" ]]; then
                echo "WARNING: empty/missing null model for ${phecode}; skipping" >&2
                printf '%s\t%s\t%s\n' "${phecode}" "-" "missing_model" >> failed_step2_singlevar.txt
                n_fail=$((n_fail + 1))
                continue
            fi
            if [[ "${have_vr}" -eq 1 && ( -z "${variance_ratio}" || ! -s "${variance_ratio}" ) ]]; then
                echo "WARNING: empty/missing variance-ratio for ${phecode}; skipping" >&2
                printf '%s\t%s\t%s\n' "${phecode}" "-" "missing_vr" >> failed_step2_singlevar.txt
                n_fail=$((n_fail + 1))
                continue
            fi
            for chromo in "${chromo_list[@]}"; do
                step2_SPAtests.R \
                    --vcfFile=~{vcf} \
                    --vcfFileIndex=~{vcf_csi} \
                    --vcfField=~{vcffield} \
                    --AlleleOrder=ref-first \
                    --chrom=${chromo} \
                    --SAIGEOutputFile=~{output_prefix}_step2Out_${phecode}_${chromo}_singlevar \
                    --minMAF=~{minimal_af} \
                    --minMAC=~{min_mac} \
                    --maxMissing=~{max_missing} \
                    --impute_method=~{impute_method} \
                    ~{"--sparseGRMFile=" + sparseGRM} \
                    ~{"--sparseGRMSampleIDFile=" + sparseGRM_IDlist} \
                    --cateVarRatioMinMACVecExclude=~{cate_var_ratio_min_mac_exclude} \
                    --cateVarRatioMaxMACVecInclude=~{cate_var_ratio_max_mac_include} \
                    --GMMATmodelFile="${model_file}" \
                    ${variance_ratio:+--varianceRatioFile="${variance_ratio}"} \
                    ~{true="--is_Firth_beta=TRUE" false="--is_Firth_beta=FALSE" is_Firth_beta} \
                    ~{true="--is_fastTest=TRUE" false="--is_fastTest=FALSE" is_fastTest} \
                    --LOCO=FALSE \
                    ~{true="--is_overrideflip=TRUE" false="" is_overridefilp} \
                    --is_output_moreDetails=TRUE
                rc=$?
                if [[ "${rc}" -ne 0 ]]; then
                    echo "WARNING: step2 singlevar failed phecode=${phecode} chrom=${chromo} exit=${rc}; continuing" >&2
                    printf '%s\t%s\t%s\n' "${phecode}" "${chromo}" "${rc}" >> failed_step2_singlevar.txt
                    n_fail=$((n_fail + 1))
                else
                    n_ok=$((n_ok + 1))
                fi
            done
        done
        echo "Step 2 singlevar kept ${n_ok} phecode×chrom runs; ${n_fail} failed (see failed_step2_singlevar.txt)" >&2
        if [[ "${n_ok}" -eq 0 ]]; then
            echo "ERROR: every step2 singlevar run failed" >&2
            exit 1
        fi
    >>>

    output {
        Array[File] outputfiles = glob("~{output_prefix}_step2Out_*_singlevar*")
        File        failed_list = "failed_step2_singlevar.txt"
    }

    runtime {
        docker:      saige_docker
        memory:      memory
        cpu:         cpu
        disks:       disk
        preemptible: preemptible
    }
}


# ---------------------------------------------------------------------------
# Task: step 2 gene-set test for every phecode in phecode_list
# ---------------------------------------------------------------------------
task RunStep2_geneset {
    meta {
        description: "Heteroplasmic variant association (all heteroplasmies) using SAIGE gene-set test — batch mode."
    }

    input {
        File vcf
        File vcf_csi
        Array[File]? variance_ratios
        Array[File] GMMATmodelFiles
        String vcffield
        Array[String] chromo_list
        Array[String] phecode_list
        String output_prefix
        Float minimal_af
        Float min_mac
        Float max_missing
        File GroupFile
        # See RunStep2_singlevariant — Step 1 fits with the sparse GRM, so Step 2 needs it too.
        File? sparseGRM
        File? sparseGRM_IDlist
        Boolean is_Firth_beta = true
        String impute_method = "mean"
        String cate_var_ratio_min_mac_exclude = "10,20.5"
        String cate_var_ratio_max_mac_include = "20.5"

        # Runtime parameters
        String memory
        String saige_docker
        Int cpu
        String disk
        Int preemptible
    }

    command <<<
        # No -e: one phecode/chrom SAIGE crash must not abort the batch.
        set -uxo pipefail

        export TMPDIR="${PWD}/tmp" TMP="${PWD}/tmp" TEMP="${PWD}/tmp"
        mkdir -p "${TMPDIR}"
        : > failed_step2_geneset.txt
        n_ok=0
        n_fail=0

        phecodes=(~{sep=' ' phecode_list})
        variance_ratios=(~{sep=' ' select_first([variance_ratios, []])})
        model_files=(~{sep=' ' GMMATmodelFiles})
        chromo_list=(~{sep=' ' chromo_list})
        have_vr=0
        if [[ ${#variance_ratios[@]} -gt 0 ]]; then
            have_vr=1
        else
            echo "WARNING: no variance-ratio files; Step 2 geneset will omit --varianceRatioFile" >&2
        fi

        # Same name-based pairing as RunStep2_singlevariant — do not zip by index.
        find_for() {
            local want="$1"; shift
            local f
            for f in "$@"; do
                if [[ "$(basename "${f}")" == *"${want}" && -s "${f}" ]]; then
                    printf '%s' "${f}"
                    return 0
                fi
            done
            return 1
        }

        if [[ ${#phecodes[@]} -eq 0 ]]; then
            echo "ERROR: phecode_list is empty" >&2
            exit 1
        fi
        if [[ ${#model_files[@]} -eq 0 ]]; then
            echo "ERROR: null_model_rdas_list is empty" >&2
            exit 1
        fi

        for i in "${!phecodes[@]}"; do
            phecode="${phecodes[$i]}"
            model_file="$(find_for "_step1Out_${phecode}.rda" "${model_files[@]}" || true)"
            variance_ratio=""
            if [[ -z "${model_file}" && ${#model_files[@]} -eq ${#phecodes[@]} ]]; then
                model_file="${model_files[$i]}"
                echo "WARNING: no .rda named *_step1Out_${phecode}.rda; falling back to position ${i} ($(basename "${model_file}")) — verify this pairing" >&2
            fi
            if [[ "${have_vr}" -eq 1 ]]; then
                variance_ratio="$(find_for "_step1Out_${phecode}.varianceRatio.txt" "${variance_ratios[@]}" || true)"
                if [[ -z "${variance_ratio}" && ${#variance_ratios[@]} -eq ${#phecodes[@]} ]]; then
                    variance_ratio="${variance_ratios[$i]}"
                    echo "WARNING: no varianceRatio named for ${phecode}; falling back to position ${i}" >&2
                fi
            fi
            if [[ -z "${model_file}" || ! -s "${model_file}" ]]; then
                echo "WARNING: empty/missing null model for ${phecode}; skipping" >&2
                printf '%s\t%s\t%s\n' "${phecode}" "-" "missing_model" >> failed_step2_geneset.txt
                n_fail=$((n_fail + 1))
                continue
            fi
            if [[ "${have_vr}" -eq 1 && ( -z "${variance_ratio}" || ! -s "${variance_ratio}" ) ]]; then
                echo "WARNING: empty/missing variance-ratio for ${phecode}; skipping" >&2
                printf '%s\t%s\t%s\n' "${phecode}" "-" "missing_vr" >> failed_step2_geneset.txt
                n_fail=$((n_fail + 1))
                continue
            fi
            for chromo in "${chromo_list[@]}"; do
                step2_SPAtests.R \
                    --vcfFile=~{vcf} \
                    --vcfFileIndex=~{vcf_csi} \
                    --vcfField=~{vcffield} \
                    --AlleleOrder=ref-first \
                    --chrom=${chromo} \
                    --SAIGEOutputFile=~{output_prefix}_step2Out_${phecode}_${chromo}_geneset \
                    --minMAF=~{minimal_af} \
                    --minMAC=~{min_mac} \
                    --maxMissing=~{max_missing} \
                    --impute_method=~{impute_method} \
                    ~{"--sparseGRMFile=" + sparseGRM} \
                    ~{"--sparseGRMSampleIDFile=" + sparseGRM_IDlist} \
                    --cateVarRatioMinMACVecExclude=~{cate_var_ratio_min_mac_exclude} \
                    --cateVarRatioMaxMACVecInclude=~{cate_var_ratio_max_mac_include} \
                    --GMMATmodelFile="${model_file}" \
                    ${variance_ratio:+--varianceRatioFile="${variance_ratio}"} \
                    ~{true="--is_Firth_beta=TRUE" false="--is_Firth_beta=FALSE" is_Firth_beta} \
                    --LOCO=FALSE \
                    --groupFile="~{GroupFile}" \
                    --annotation_in_groupTest="non_coding_transcript_exon,start_lost,stop_gained;start_lost,missense,synonymous,frameshift,stop_gained,stop_lost,stop_retained,dloop" \
                    --maxMAF_in_groupTest=0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5 \
                    --lr_PCutoffforFirth=0.05 \
                    --is_output_markerList_in_groupTest=TRUE \
                    --is_output_moreDetails=TRUE
                rc=$?
                if [[ "${rc}" -ne 0 ]]; then
                    echo "WARNING: step2 geneset failed phecode=${phecode} chrom=${chromo} exit=${rc}; continuing" >&2
                    printf '%s\t%s\t%s\n' "${phecode}" "${chromo}" "${rc}" >> failed_step2_geneset.txt
                    n_fail=$((n_fail + 1))
                else
                    n_ok=$((n_ok + 1))
                fi
            done
        done
        echo "Step 2 geneset kept ${n_ok} phecode×chrom runs; ${n_fail} failed (see failed_step2_geneset.txt)" >&2
        if [[ "${n_ok}" -eq 0 ]]; then
            echo "ERROR: every step2 geneset run failed" >&2
            exit 1
        fi
    >>>

    output {
        Array[File] outputfiles = glob("~{output_prefix}_step2Out_*_geneset*")
        File        failed_list = "failed_step2_geneset.txt"
    }

    runtime {
        docker:      saige_docker
        memory:      memory
        cpu:         cpu
        disks:       disk
        preemptible: preemptible
    }
}
