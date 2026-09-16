version 1.0

workflow SAIGE_phewas {
    meta {
        description: "PheWAS workflow for variant association using SAIGE — each task processes a batch of phecodes."
    }

    input {
        File plink_bed_file
        File plink_bim_file
        File plink_fam_file
        File sparseGRM
        File sparseGRM_IDlist
        Array[String] phecode_list
        File homoplasmic_vcf
        File homoplasmic_vcf_csi
        File? heteroplasmic_vcf
        File? heteroplasmic_vcf_csi
        File? GroupFile
        String chromosome
        File phenotype_file
        String trait_type
        String output_prefix
        String saige_docker
        String vcffield
        Boolean is_overridefilp

        String memory = "8G"
        
        Float single_variant_minimal_af = 0.01
        Int single_variant_min_mac = 20
        Float gene_set_minimal_af = 0
        Float gene_set_min_mac = 0.5

        # MAC category bounds for categorical variance ratios, shared by Step 1 and Step 2 so the
        # two cannot drift. SAIGE requires length(min_exclude) == length(max_include) + 1, and
        # length(min_exclude) == the number of "null" rows Step 1 writes to the varianceRatio file.
        # If these disagree between steps, SAIGE applies each ratio to the wrong MAC bin silently.
        String cate_var_ratio_min_mac_exclude = "100,500"
        String cate_var_ratio_max_mac_include = "500"

        Int cpu = 2
        # 50G is tight: AoU sparse GRM + plink + 15 null .rda files. Empty Step 1
        # outputs are often a full /tmp or a full local-disk, not a SAIGE bug.
        String disk_size = "local-disk 200 HDD"
        Int preemptible = 1
    }
    
    # Step 1: fit one null GLMM per phecode in a single batch task
    call RunFitNullGLMM {
        input:
            plink_bed_file    = plink_bed_file,
            plink_bim_file    = plink_bim_file,
            plink_fam_file    = plink_fam_file,
            sparseGRM         = sparseGRM,
            sparseGRM_IDlist  = sparseGRM_IDlist,
            phenotype_file    = phenotype_file,
            phecode_list      = phecode_list,
            trait_type        = trait_type,
            output_prefix     = output_prefix,
            saige_docker      = saige_docker,
            cate_var_ratio_min_mac_exclude = cate_var_ratio_min_mac_exclude,
            cate_var_ratio_max_mac_include = cate_var_ratio_max_mac_include,
            memory            = memory,
            cpu               = cpu,
            disk              = disk_size,
            preemptible       = preemptible
    }

    # Step 2a: single-variant test (homoplasmic VCF) for all phecodes
    call RunStep2_singlevariant {
        input:
            vcf              = homoplasmic_vcf,
            vcf_csi          = homoplasmic_vcf_csi,
            variance_ratios  = RunFitNullGLMM.variance_ratio_txts,
            GMMATmodelFiles  = RunFitNullGLMM.null_model_rdas,
            sparseGRM        = sparseGRM,
            sparseGRM_IDlist = sparseGRM_IDlist,
            vcffield         = vcffield,
            chromo           = chromosome,
            phecode_list     = read_lines(RunFitNullGLMM.successful_phecode_list),
            output_prefix    = output_prefix,
            minimal_af       = single_variant_minimal_af,
            min_mac          = single_variant_min_mac,
            is_overridefilp  = is_overridefilp,
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
                variance_ratios  = RunFitNullGLMM.variance_ratio_txts,
                GMMATmodelFiles  = RunFitNullGLMM.null_model_rdas,
                sparseGRM        = sparseGRM,
                sparseGRM_IDlist = sparseGRM_IDlist,
                vcffield         = vcffield,
                chromo           = chromosome,
                phecode_list     = read_lines(RunFitNullGLMM.successful_phecode_list),
                output_prefix    = output_prefix,
                minimal_af       = gene_set_minimal_af,
                min_mac          = gene_set_min_mac,
                GroupFile        = select_first([GroupFile, ""]),
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
        Array[File] null_model_rdas        = RunFitNullGLMM.null_model_rdas
        Array[File] variance_ratio_txts    = RunFitNullGLMM.variance_ratio_txts
        Array[File] singlevariant_outputs  = RunStep2_singlevariant.outputfiles
        Array[File]? geneset_outputs        = RunStep2_geneset.outputfiles
    }
}


# ---------------------------------------------------------------------------
# Task: fit null GLMM for every phecode in phecode_list
# ---------------------------------------------------------------------------
task RunFitNullGLMM {

    input {
        File plink_bed_file
        File plink_bim_file
        File plink_fam_file
        File sparseGRM
        File sparseGRM_IDlist
        File phenotype_file
        Array[String] phecode_list
        Array[String] covariate_list = ["PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","sex","age","age2","age_sex","age2_sex","coverage", "GC"]
        Array[String] categorical_covariate = ["sex", "GC"]
        String trait_type
        String output_prefix
        String saige_docker

        # Must match the values passed to step2_SPAtests.R. See the workflow-level inputs.
        String cate_var_ratio_min_mac_exclude = "10,20.5"
        String cate_var_ratio_max_mac_include = "20.5"

        # Runtime parameters
        String memory
        Int cpu
        String disk
        Int preemptible
    }

    String inv_norm_arg = if (trait_type == "quantitative") then "--invNormalize=TRUE" else ""

    command <<<
        set -uxo pipefail

        # R tempfile()/save() default to /tmp, which on Terra/dataproc is a small
        # boot disk. A failed write leaves a 0-byte .rda that glob() still collects.
        export TMPDIR="${PWD}/tmp" TMP="${PWD}/tmp" TEMP="${PWD}/tmp"
        mkdir -p "${TMPDIR}"
        df -h . /tmp >&2 || true

        phecodes=(~{sep=' ' phecode_list})
        successful_phecodes=()
        declare -A seen=()
        unique_phecodes=()

        # Cromwell glob() is unordered and unique-by-basename, so a duplicate phecode
        # overwrites the same .rda/.varianceRatio.txt but is recorded twice in the
        # success list — Step 2 then has more names than files.
        for phecode in "${phecodes[@]}"; do
            if [[ -n "${seen[$phecode]+x}" ]]; then
                echo "WARNING: duplicate phecode ${phecode} in phecode_list, skipping" >&2
                continue
            fi
            seen[$phecode]=1
            unique_phecodes+=("${phecode}")
        done

        for phecode in "${unique_phecodes[@]}"; do
            model_file="~{output_prefix}_step1Out_${phecode}.rda"
            vr_file="~{output_prefix}_step1Out_${phecode}.varianceRatio.txt"
            # SAIGE often writes .rda + .varianceRatio.txt and then exits non-zero
            # (R warnings / glm.fit). Success is the files, not the R exit code.
            step1_fitNULLGLMM.R \
                --bedFile="~{plink_bed_file}" \
                --bimFile="~{plink_bim_file}" \
                --famFile="~{plink_fam_file}" \
                --useSparseGRMtoFitNULL=TRUE \
                --useSparseGRMforVarRatio=TRUE \
                --sparseGRMFile="~{sparseGRM}" \
                --sparseGRMSampleIDFile="~{sparseGRM_IDlist}" \
                --phenoFile="~{phenotype_file}" \
                --phenoCol="${phecode}" \
                --covarColList="~{sep=',' covariate_list}" \
                --qCovarColList="~{sep=',' categorical_covariate}" \
                --sampleIDColinphenoFile=person_id \
                --traitType=~{trait_type} \
                --isCateVarianceRatio=TRUE \
                --cateVarRatioMinMACVecExclude=~{cate_var_ratio_min_mac_exclude} \
                --cateVarRatioMaxMACVecInclude=~{cate_var_ratio_max_mac_include} \
                --IsOverwriteVarianceRatioFile=TRUE \
                --outputPrefix="~{output_prefix}_step1Out_${phecode}" \
                ~{inv_norm_arg}
            rc=$?
            # -s: exist AND non-empty. SAIGE creates 0-byte placeholders, then
            # save() can fail (full /tmp, full local-disk) and leave them empty.
            if [[ -s "${model_file}" && -s "${vr_file}" ]]; then
                if [[ "${rc}" -ne 0 ]]; then
                    echo "WARNING: phecode ${phecode} step1 exited ${rc} but outputs exist; keeping" >&2
                fi
                successful_phecodes+=("${phecode}")
            else
                echo "WARNING: phecode ${phecode} failed (exit ${rc}); missing or empty ${model_file} ($(stat -c%s "${model_file}" 2>/dev/null || echo 0) bytes) / ${vr_file} ($(stat -c%s "${vr_file}" 2>/dev/null || echo 0) bytes)" >&2
            fi
        done
        echo "Step 1 kept ${#successful_phecodes[@]} / ${#unique_phecodes[@]} phecodes" >&2

        if [[ ${#successful_phecodes[@]} -gt 0 ]]; then
            printf '%s\n' "${successful_phecodes[@]}" > successful_phecodes.txt
        else
            : > successful_phecodes.txt
            echo "WARNING: no phecodes produced Step 1 outputs" >&2
        fi
    >>>

    output {
        File        successful_phecode_list = "successful_phecodes.txt"
        Array[File] null_model_rdas         = glob("~{output_prefix}_step1Out_*.rda")
        Array[File] variance_ratio_txts     = glob("~{output_prefix}_step1Out_*.varianceRatio.txt")
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
# Task: step 2 single-variant SPA test for every phecode in phecode_list
# ---------------------------------------------------------------------------
task RunStep2_singlevariant {

    input {
        File vcf
        File vcf_csi
        Array[File] variance_ratios
        Array[File] GMMATmodelFiles
        # Step 1 in this workflow always fits with --useSparseGRMtoFitNULL=TRUE, so Step 2 must be
        # given the same sparse GRM — otherwise SAIGE silently falls back to the "null" variance
        # ratios instead of the sparse-GRM ones (the guard for this is commented out upstream).
        File? sparseGRM
        File? sparseGRM_IDlist
        String vcffield
        String chromo
        Array[String] phecode_list
        String output_prefix
        Boolean is_overridefilp = false
        Boolean is_fastTest = true
        Float minimal_af
        Int min_mac
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
        set -euxo pipefail

        export TMPDIR="${PWD}/tmp" TMP="${PWD}/tmp" TEMP="${PWD}/tmp"
        mkdir -p "${TMPDIR}"

        phecodes=(~{sep=' ' phecode_list})
        variance_ratios=(~{sep=' ' variance_ratios})
        model_files=(~{sep=' ' GMMATmodelFiles})

        # glob() order is not the success-list order. Match by exact basename so
        # e.g. EM_200 does not pick up EM_200.1, and a length mismatch cannot
        # silently pair the wrong null model.
        # If successful_phecodes.txt was empty (SAIGE exited non-zero after
        # writing files), recover names from the localized .rda files.
        expected_prefix="~{output_prefix}_step1Out_"
        if [[ ${#phecodes[@]} -eq 0 ]]; then
            echo "WARNING: successful_phecode_list was empty; inferring phecodes from Step 1 .rda files" >&2
            for f in "${model_files[@]+"${model_files[@]}"}"; do
                [[ -s "${f}" ]] || continue
                stem="$(basename "${f}" .rda)"
                if [[ "${stem}" == "${expected_prefix}"* ]]; then
                    phecodes+=("${stem#"${expected_prefix}"}")
                fi
            done
        fi
        if [[ ${#phecodes[@]} -eq 0 ]]; then
            echo "ERROR: no phecodes to test (success list empty and no inferable .rda files)" >&2
            exit 1
        fi
        if [[ ${#model_files[@]} -eq 0 || ${#variance_ratios[@]} -eq 0 ]]; then
            echo "ERROR: Step 1 produced no model or variance-ratio files" >&2
            exit 1
        fi

        for phecode in "${phecodes[@]}"; do
            expected_rda="~{output_prefix}_step1Out_${phecode}.rda"
            expected_vr="~{output_prefix}_step1Out_${phecode}.varianceRatio.txt"
            model_file=""
            variance_ratio=""
            for f in "${model_files[@]}"; do
                if [[ "$(basename "${f}")" == "${expected_rda}" && -s "${f}" ]]; then
                    model_file="${f}"
                    break
                fi
            done
            for f in "${variance_ratios[@]}"; do
                if [[ "$(basename "${f}")" == "${expected_vr}" && -s "${f}" ]]; then
                    variance_ratio="${f}"
                    break
                fi
            done
            if [[ -z "${model_file}" || -z "${variance_ratio}" ]]; then
                echo "ERROR: no Step 1 outputs matching ${phecode}" >&2
                echo "  expected rda: ${expected_rda}" >&2
                echo "  expected vr:  ${expected_vr}" >&2
                exit 1
            fi

            step2_SPAtests.R \
                --vcfFile=~{vcf} \
                --vcfFileIndex=~{vcf_csi} \
                --vcfField=~{vcffield} \
                --AlleleOrder=ref-first \
                --chrom=~{chromo} \
                --SAIGEOutputFile=~{output_prefix}_step2Out_${phecode}_singlevar \
                --minMAF=~{minimal_af} \
                --minMAC=~{min_mac} \
                ~{"--sparseGRMFile=" + sparseGRM} \
                ~{"--sparseGRMSampleIDFile=" + sparseGRM_IDlist} \
                --cateVarRatioMinMACVecExclude=~{cate_var_ratio_min_mac_exclude} \
                --cateVarRatioMaxMACVecInclude=~{cate_var_ratio_max_mac_include} \
                --GMMATmodelFile="${model_file}" \
                --varianceRatioFile="${variance_ratio}" \
                --is_Firth_beta=TRUE \
                ~{true="--is_fastTest=TRUE" false="--is_fastTest=FALSE" is_fastTest} \
                --LOCO=FALSE \
                ~{true="--is_overrideflip=TRUE" false="" is_overridefilp} \
                --is_output_moreDetails=TRUE
        done
    >>>

    output {
        Array[File] outputfiles = glob("~{output_prefix}_step2Out_*_singlevar*")
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
        Array[File] variance_ratios
        Array[File] GMMATmodelFiles
        String vcffield
        String chromo
        Array[String] phecode_list
        String output_prefix
        Float minimal_af
        Float min_mac
        File GroupFile
        # See RunStep2_singlevariant — Step 1 fits with the sparse GRM, so Step 2 needs it too.
        File? sparseGRM
        File? sparseGRM_IDlist
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
        set -euxo pipefail

        export TMPDIR="${PWD}/tmp" TMP="${PWD}/tmp" TEMP="${PWD}/tmp"
        mkdir -p "${TMPDIR}"

        phecodes=(~{sep=' ' phecode_list})
        variance_ratios=(~{sep=' ' variance_ratios})
        model_files=(~{sep=' ' GMMATmodelFiles})

        # Same basename lookup as RunStep2_singlevariant — do not zip by index.
        expected_prefix="~{output_prefix}_step1Out_"
        if [[ ${#phecodes[@]} -eq 0 ]]; then
            echo "WARNING: successful_phecode_list was empty; inferring phecodes from Step 1 .rda files" >&2
            for f in "${model_files[@]+"${model_files[@]}"}"; do
                [[ -s "${f}" ]] || continue
                stem="$(basename "${f}" .rda)"
                if [[ "${stem}" == "${expected_prefix}"* ]]; then
                    phecodes+=("${stem#"${expected_prefix}"}")
                fi
            done
        fi
        if [[ ${#phecodes[@]} -eq 0 ]]; then
            echo "ERROR: no phecodes to test (success list empty and no inferable .rda files)" >&2
            exit 1
        fi
        if [[ ${#model_files[@]} -eq 0 || ${#variance_ratios[@]} -eq 0 ]]; then
            echo "ERROR: Step 1 produced no model or variance-ratio files" >&2
            exit 1
        fi

        for phecode in "${phecodes[@]}"; do
            expected_rda="~{output_prefix}_step1Out_${phecode}.rda"
            expected_vr="~{output_prefix}_step1Out_${phecode}.varianceRatio.txt"
            model_file=""
            variance_ratio=""
            for f in "${model_files[@]}"; do
                if [[ "$(basename "${f}")" == "${expected_rda}" && -s "${f}" ]]; then
                    model_file="${f}"
                    break
                fi
            done
            for f in "${variance_ratios[@]}"; do
                if [[ "$(basename "${f}")" == "${expected_vr}" && -s "${f}" ]]; then
                    variance_ratio="${f}"
                    break
                fi
            done
            if [[ -z "${model_file}" || -z "${variance_ratio}" ]]; then
                echo "ERROR: no Step 1 outputs matching ${phecode}" >&2
                echo "  expected rda: ${expected_rda}" >&2
                echo "  expected vr:  ${expected_vr}" >&2
                exit 1
            fi

            step2_SPAtests.R \
                --vcfFile=~{vcf} \
                --vcfFileIndex=~{vcf_csi} \
                --vcfField=~{vcffield} \
                --AlleleOrder=ref-first \
                --chrom=~{chromo} \
                --SAIGEOutputFile=~{output_prefix}_step2Out_${phecode}_geneset \
                --minMAF=~{minimal_af} \
                --minMAC=~{min_mac} \
                ~{"--sparseGRMFile=" + sparseGRM} \
                ~{"--sparseGRMSampleIDFile=" + sparseGRM_IDlist} \
                --cateVarRatioMinMACVecExclude=~{cate_var_ratio_min_mac_exclude} \
                --cateVarRatioMaxMACVecInclude=~{cate_var_ratio_max_mac_include} \
                --GMMATmodelFile="${model_file}" \
                --varianceRatioFile="${variance_ratio}" \
                --is_Firth_beta=TRUE \
                --LOCO=FALSE \
                --groupFile="~{GroupFile}" \
                --annotation_in_groupTest="non_coding_transcript_exon,start_lost,stop_gained;start_lost,missense,synonymous,frameshift,stop_gained,stop_lost,stop_retained,dloop" \
                --maxMAF_in_groupTest=0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5 \
                --lr_PCutoffforFirth=0.05 \
                --is_output_markerList_in_groupTest=TRUE \
                --is_output_moreDetails=TRUE
        done
    >>>

    output {
        Array[File] outputfiles = glob("~{output_prefix}_step2Out_*_geneset*")
    }

    runtime {
        docker:      saige_docker
        memory:      memory
        cpu:         cpu
        disks:       disk
        preemptible: preemptible
    }
}
