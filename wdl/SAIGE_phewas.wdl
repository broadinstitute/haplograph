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
        Array[File]? null_model_rdas_list
        Array[File]? variance_ratio_txts_list
        Array[String] chromosome_list
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

        # When false, Step 2 omits --varianceRatioFile and sets --is_fastTest=FALSE
        # (SAIGE uses variance ratio = 1). Step 1 still writes the VR files.
        Boolean use_varianceRatio = true

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

        # A phecode is treated as sex-specific when at least this fraction of its cases are one
        # sex. Not 1.0: AoU sex_at_birth disagrees with the clinical code for a handful of
        # participants, and under an exact rule a single discordant case leaves the phenotype
        # unrestricted — exactly the inflation this guard exists to prevent (GU_619.2 ovarian
        # cyst). Discordant cases are dropped with the rest of the other sex.
        Float sex_specific_case_fraction = 0.99

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
    
    # Step 1: fit one null GLMM per phecode in a single batch task
    if (!defined(null_model_rdas_list)) { 
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
                covariate_list = covariate_list,
                categorical_covariate = categorical_covariate_list,
                sex_specific_case_fraction = sex_specific_case_fraction,

                memory            = step1_memory,
                cpu               = cpu,
                disk              = disk_size,
                preemptible       = preemptible
        }
    }
    # When Step 1 ran, use its successful subset; otherwise use the workflow phecode_list.
    Array[String] step2_phecode_list = if defined(RunFitNullGLMM.successful_phecode_list)
        then read_lines(select_first([RunFitNullGLMM.successful_phecode_list]))
        else phecode_list

    # Null models are always required (user list or Step 1). VR files are optional.
    # use_varianceRatio=false still localizes any VR files that exist, but Step 2
    # ignores them: is_fastTest=FALSE and no --varianceRatioFile (SAIGE uses ratio=1).
    # The same happens when Step 1 is skipped and variance_ratio_txts_list is omitted.
    Array[File] step2_null_models = select_first([null_model_rdas_list, RunFitNullGLMM.null_model_rdas])
    Array[File]? step2_variance_ratios = if defined(variance_ratio_txts_list)
        then variance_ratio_txts_list
        else RunFitNullGLMM.variance_ratio_txts

    # Step 2a: single-variant test (homoplasmic VCF) for all phecodes
    call RunStep2_singlevariant {
        input:
            vcf              = homoplasmic_vcf,
            vcf_csi          = homoplasmic_vcf_csi,
            variance_ratios  = step2_variance_ratios,
            GMMATmodelFiles  = step2_null_models,
            sparseGRM        = sparseGRM,
            sparseGRM_IDlist = sparseGRM_IDlist,
            vcffield         = vcffield,
            chromo_list           = chromosome_list,
            phecode_list     = step2_phecode_list,
            output_prefix    = output_prefix,
            minimal_af       = single_variant_minimal_af,
            min_mac          = single_variant_min_mac,
            max_missing      = single_variant_max_missing,
            is_overridefilp  = is_overridefilp,
            is_Firth_beta    = is_Firth_beta,
            use_varianceRatio = use_varianceRatio,
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
                variance_ratios  = step2_variance_ratios,
                GMMATmodelFiles  = step2_null_models,
                sparseGRM        = sparseGRM,
                sparseGRM_IDlist = sparseGRM_IDlist,
                vcffield         = vcffield,
                chromo_list           = chromosome_list,
                phecode_list     = step2_phecode_list,
                output_prefix    = output_prefix,
                minimal_af       = gene_set_minimal_af,
                min_mac          = gene_set_min_mac,
                max_missing      = gene_set_max_missing,
                GroupFile        = select_first([GroupFile, ""]),
                is_Firth_beta    = is_Firth_beta,
                use_varianceRatio = use_varianceRatio,
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
        File?       covariate_plan          = RunFitNullGLMM.covariate_plan
        Array[File]? null_model_rdas        = RunFitNullGLMM.null_model_rdas
        Array[File]? variance_ratio_txts    = RunFitNullGLMM.variance_ratio_txts
        Array[File] singlevariant_outputs  = RunStep2_singlevariant.outputfiles
        File        singlevariant_failed   = RunStep2_singlevariant.failed_list
        Array[File]? geneset_outputs        = RunStep2_geneset.outputfiles
        File?       geneset_failed          = RunStep2_geneset.failed_list
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
        Array[String] covariate_list
        Array[String] categorical_covariate
        String trait_type
        String output_prefix
        String saige_docker

        # See the workflow-level input of the same name.
        Float sex_specific_case_fraction = 0.99

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

        # A sex-specific phecode (every case is one sex; the other sex is a control)
        # is fit on the case sex only. Leaving the other sex in makes case status a
        # stand-in for sex, and the score test inflates even with sex in the model.
        # After that restriction, sex is constant and age_sex / age2_sex are copies
        # of age / age2, so those columns are dropped. A covariate that is constant
        # or an exact copy on the samples actually fit is dropped the same way
        # (a rare genome-center dummy can go all-zero). The covariate list is shared
        # by every phenotype, so this has to happen per phecode.
        printf '%s\n' "${unique_phecodes[@]}" > requested_phecodes.txt
        cat > drop_constant_covars.R << 'ENDRS'
args <- commandArgs(trailingOnly = TRUE)
pheno_path <- args[[1]]
phecode_path <- args[[2]]
covars <- strsplit(args[[3]], ",", fixed = TRUE)[[1]]
covars <- covars[nzchar(covars)]
qcovars <- strsplit(args[[4]], ",", fixed = TRUE)[[1]]
qcovars <- qcovars[nzchar(qcovars)]
out_path <- args[[5]]
trait_type <- if (length(args) >= 6) args[[6]] else "binary"
sex_case_frac <- if (length(args) >= 7) suppressWarnings(as.numeric(args[[7]])) else 0.99
if (!is.finite(sex_case_frac) || sex_case_frac <= 0 || sex_case_frac > 1) sex_case_frac <- 0.99
if (length(covars) < 1) stop("No covariates listed")

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("data.table is required to screen covariates")
}
phecodes <- readLines(phecode_path)
phecodes <- phecodes[nzchar(phecodes)]
hdr <- names(data.table::fread(pheno_path, nrows = 0))
id_col <- "person_id"
if (!(id_col %in% hdr)) stop("Phenotype file has no person_id column")
missing_cov <- setdiff(covars, hdr)
if (length(missing_cov)) {
  stop("Phenotype file is missing covariates: ", paste(missing_cov, collapse = ","))
}
present <- phecodes[phecodes %in% hdr]
# person_id must stay text. AoU IDs pasted through numeric lose digits, and
# SampleIDIncludeFile is matched as character against the phenotype file.
dt <- data.table::fread(
  pheno_path,
  select = unique(c(id_col, present, covars)),
  colClasses = list(character = id_col)
)
# A covariate named in qCovarColList, or one whose values are labels rather than
# numbers (GC = "BI"), is categorical. SAIGE dummy-codes those via qCovarColList.
# A numeric column with a few bad tokens stays numeric; those tokens become NA.
is_cat <- setNames(covars %in% qcovars, covars)
for (cv in covars) {
  raw <- trimws(as.character(dt[[cv]]))
  raw[is.na(raw) | !nzchar(raw) | raw %in% c("NA", "NaN", "nan", ".")] <- NA
  num <- suppressWarnings(as.numeric(raw))
  present <- !is.na(raw)
  bad <- present & is.na(num)
  if (isTRUE(is_cat[[cv]]) || (any(present) && mean(bad[present]) > 0.5)) {
    is_cat[[cv]] <- TRUE
    if (!(cv %in% qcovars)) {
      qcovars <- c(qcovars, cv)
      message(sprintf("Treating %s as categorical (example: %s)", cv, raw[which(bad)[1]]))
    }
    dt[[cv]] <- raw
  } else {
    if (any(bad)) {
      message(sprintf(
        "Setting %d non-numeric %s value(s) to NA (example: %s)",
        sum(bad), cv, raw[which(bad)[1]]
      ))
    }
    dt[[cv]] <- num
  }
}

# Keep a covariate unless it is constant, or an exact copy of one already kept
# (sex fixed at 1 makes age_sex identical to age). Categorical columns are
# constant when they have a single level; they are not checked for correlation.
screen_covars <- function(sub, covars, is_cat) {
  constant <- character(0)
  numeric_cov <- covars[!is_cat[covars]]
  cat_cov <- covars[is_cat[covars]]
  if (length(numeric_cov)) {
    sds <- vapply(numeric_cov, function(cv) stats::sd(sub[[cv]]), numeric(1))
    constant <- c(constant, numeric_cov[!is.finite(sds) | sds == 0])
  }
  if (length(cat_cov)) {
    nlev <- vapply(cat_cov, function(cv) {
      v <- as.character(sub[[cv]])
      v <- v[!is.na(v) & nzchar(v)]
      length(unique(v))
    }, integer(1))
    constant <- c(constant, cat_cov[nlev < 2])
  }
  remain <- setdiff(covars, constant)
  aliased <- character(0)
  kept <- character(0)
  for (cv in remain) {
    if (isTRUE(is_cat[[cv]])) {
      kept <- c(kept, cv)
      next
    }
    x <- sub[[cv]]
    is_alias <- FALSE
    for (prev in kept) {
      if (isTRUE(is_cat[[prev]])) next
      r <- suppressWarnings(stats::cor(x, sub[[prev]]))
      if (is.finite(r) && abs(r) > 1 - 1e-8) {
        is_alias <- TRUE
        break
      }
    }
    if (is_alias) aliased <- c(aliased, cv) else kept <- c(kept, cv)
  }
  list(keep = kept, dropped = c(constant, aliased))
}

# TRUE if the value is one of the tokens used for "no value".
is_missing_token <- function(v) {
  v <- trimws(as.character(v))
  is.na(v) | !nzchar(v) | v %in% c("NA", "NaN", "nan", ".")
}

# Phecode values as 0/1, NA for anything else. data.table reads an all-TRUE /
# FALSE column as logical and some phenotype exports use Y/N, so a phenotype
# that is binary in substance can arrive as text; under the old numeric-only
# rule those phecodes looked non-binary and silently skipped sex restriction.
as01 <- function(y) {
  v <- trimws(as.character(y))
  v[is_missing_token(v)] <- NA
  out <- suppressWarnings(as.numeric(v))
  lv <- toupper(v)
  out[is.na(out) & !is.na(lv) & lv %in% c("TRUE", "T", "YES", "Y", "CASE")] <- 1
  out[is.na(out) & !is.na(lv) & lv %in% c("FALSE", "F", "NO", "N", "CONTROL")] <- 0
  out[!is.na(out) & !(abs(out) < 1e-8 | abs(out - 1) < 1e-8)] <- NA
  out
}

# Binary phecode whose cases are (almost) all one sex, with the other sex
# present as controls. Returns that case-sex value, or NA when the phecode is
# both-sex, already single-sex, or not coded 0/1.
#
# `frac` is deliberately below 1. Requiring every case to share a sex means one
# discordant row — AoU sex_at_birth vs the clinical code, a data-entry error —
# leaves a female-only phenotype fit on both sexes, where case status stands in
# for sex and any sex-correlated marker (chrX/chrY, or a sex-differential
# dosage call) separates cases from controls outright. That is the mechanism
# behind lambda in the hundreds, not a few percent of inflation.
sex_to_keep <- function(y01, sex, frac = 0.99) {
  is_case <- !is.na(y01) & abs(y01 - 1) < 1e-8
  is_ctrl <- !is.na(y01) & abs(y01) < 1e-8
  if (!any(is_case)) return(NA_character_)
  known <- !is.na(sex) & nzchar(sex)
  case_sex_tab <- table(sex[is_case & known])
  if (length(case_sex_tab) < 1) return(NA_character_)
  top <- names(case_sex_tab)[which.max(case_sex_tab)]
  if (max(case_sex_tab) / sum(case_sex_tab) < frac) return(NA_character_)
  # Only restrict when the other sex is actually in the analysis set; a
  # phenotype already limited to one sex needs no exclusion.
  other <- (is_case | is_ctrl) & known & sex != top
  if (!any(other)) return(NA_character_)
  top
}

out <- file(out_path, "w")
on.exit(close(out), add = TRUE)
emit <- function(...) cat(paste(c(...), collapse = "\t"), "\n", sep = "", file = out)
for (phe in phecodes) {
  if (!(phe %in% names(dt))) {
    emit(phe, "SKIP", "", "", "MISSING_COLUMN", "", "")
    next
  }
  y <- dt[[phe]]
  keep_rows <- !is_missing_token(y)
  include_file <- ""
  pheno_file <- ""
  # Quantitative traits are never 0/1, so sex restriction does not apply; the
  # decision is made from trait_type rather than inferred from the values.
  if (trait_type == "binary" && "sex" %in% names(dt)) {
    sex <- trimws(as.character(dt[["sex"]]))
    sex[is_missing_token(sex)] <- NA
    y01 <- as01(y)
    # A value SAIGE cannot read as 0/1 is not a control; leave such a phecode
    # alone rather than guessing at its coding.
    unreadable <- sum(keep_rows & is.na(y01))
    keep_sex <- if (unreadable > 0) NA_character_ else sex_to_keep(y01, sex, sex_case_frac)
    if (unreadable > 0) {
      message(sprintf(
        "%s: %d value(s) are neither 0 nor 1; not screening for sex-specificity",
        phe, unreadable
      ))
    }
    if (!is.na(keep_sex)) {
      in_sex <- !is.na(sex) & sex == keep_sex
      n_case <- sum(!is.na(y01) & abs(y01 - 1) < 1e-8)
      n_case_kept <- sum(!is.na(y01) & abs(y01 - 1) < 1e-8 & in_sex)
      n_excl <- sum(keep_rows & !in_sex)
      keep_rows <- keep_rows & in_sex
      ids <- unique(as.character(dt[[id_col]][in_sex]))
      ids <- ids[!is.na(ids) & nzchar(ids)]
      safe <- gsub("[^A-Za-z0-9._-]", "_", phe)
      include_file <- sprintf("include_%s.txt", safe)
      writeLines(ids, include_file)
      # Belt and braces: also hand SAIGE a phenotype file where the excluded
      # sex (and any case discordant with it) is NA, so the restriction holds
      # even on a SAIGE build without --SampleIDIncludeFile, where that flag
      # would otherwise be ignored or fatal.
      pheno_file <- sprintf("pheno_%s.tsv", safe)
      # copy(): column subsetting shallow-copies, and set() would otherwise
      # overwrite dt's own phecode column for every later phecode.
      sub_dt <- data.table::copy(dt[, unique(c(id_col, phe, covars)), with = FALSE])
      data.table::set(sub_dt, i = which(!in_sex), j = phe, value = NA)
      data.table::fwrite(sub_dt, pheno_file, sep = "\t", na = "NA", quote = FALSE)
      message(sprintf(
        paste0("%s: %d/%d cases are sex=%s; fitting on that sex only ",
               "(excluding %d samples, of which %d discordant cases)"),
        phe, n_case_kept, n_case, keep_sex, n_excl, n_case - n_case_kept
      ))
    }
  }
  sub <- dt[keep_rows, covars, with = FALSE]
  ok <- stats::complete.cases(sub)
  sub <- sub[ok, , drop = FALSE]
  if (nrow(sub) < 2) {
    emit(phe, "SKIP", "", "", "TOO_FEW_SAMPLES", "", "")
    next
  }
  scr <- screen_covars(sub, covars, is_cat)
  if (length(scr$keep) < 1) {
    emit(phe, "SKIP", "", "", "ALL_CONSTANT", "", "")
    next
  }
  keep_q <- intersect(qcovars, scr$keep)
  emit(phe, "OK",
       paste(scr$keep, collapse = ","),
       paste(keep_q, collapse = ","),
       paste(scr$dropped, collapse = ","),
       include_file,
       pheno_file)
}
ENDRS
        # saige:1.5.1.noflip keeps R inside pixi. Rscript is not on PATH;
        # step1_fitNULLGLMM.R reaches it with this same pixi invocation.
        if command -v Rscript >/dev/null 2>&1; then
            rscript=(Rscript)
        elif command -v pixi >/dev/null 2>&1 && [[ -f /app/pixi.toml ]]; then
            rscript=(pixi run --manifest-path /app/pixi.toml Rscript)
        else
            echo "ERROR: Rscript not on PATH and pixi manifest /app/pixi.toml not found" >&2
            exit 1
        fi
        if ! "${rscript[@]}" --vanilla drop_constant_covars.R \
            "~{phenotype_file}" \
            requested_phecodes.txt \
            "~{sep=',' covariate_list}" \
            "~{sep=',' categorical_covariate}" \
            covariates_by_phecode.tsv \
            "~{trait_type}" \
            "~{sex_specific_case_fraction}"
        then
            echo "ERROR: failed to screen covariates for constant columns" >&2
            exit 1
        fi

        declare -A covar_for=()
        declare -A qcovar_for=()
        declare -A dropped_for=()
        declare -A status_for=()
        declare -A include_for=()
        declare -A pheno_for=()
        while IFS=$'\t' read -r phe status covs qcovs dropped include pheno; do
            [[ -n "${phe}" ]] || continue
            status_for["${phe}"]="${status}"
            covar_for["${phe}"]="${covs}"
            qcovar_for["${phe}"]="${qcovs}"
            dropped_for["${phe}"]="${dropped}"
            include_for["${phe}"]="${include}"
            pheno_for["${phe}"]="${pheno}"
        done < covariates_by_phecode.tsv

        # --SampleIDIncludeFile exists in SAIGE 1.x step 1, but an older or
        # patched build that lacks it would either die on the unknown flag or
        # ignore it and fit both sexes anyway. Check once; the per-phecode
        # phenotype file (other sex set to NA) carries the restriction either way.
        have_include_opt=1
        if ! step1_fitNULLGLMM.R --help 2>&1 | grep -q -- "--SampleIDIncludeFile"; then
            have_include_opt=0
            echo "WARNING: this SAIGE build has no --SampleIDIncludeFile; relying on the per-phecode phenotype file to restrict sex" >&2
        fi

        for phecode in "${unique_phecodes[@]}"; do
            model_file="~{output_prefix}_step1Out_${phecode}.rda"
            vr_file="~{output_prefix}_step1Out_${phecode}.varianceRatio.txt"
            if [[ "${status_for[${phecode}]:-}" != "OK" ]]; then
                echo "WARNING: phecode ${phecode} skipped before step 1 (${dropped_for[${phecode}]:-no covariate plan})" >&2
                continue
            fi
            if [[ -n "${dropped_for[${phecode}]}" ]]; then
                echo "Dropping covariate(s) constant or perfectly collinear on ${phecode} samples: ${dropped_for[${phecode}]}" >&2
            fi
            pheno_file="~{phenotype_file}"
            if [[ -n "${pheno_for[${phecode}]:-}" && -s "${pheno_for[${phecode}]}" ]]; then
                pheno_file="${pheno_for[${phecode}]}"
                echo "Fitting ${phecode} on the sex that has cases (${pheno_file}; the other sex is NA)" >&2
            fi
            # SAIGE often writes .rda + .varianceRatio.txt and then exits non-zero
            # (R warnings / glm.fit). Success is the files, not the R exit code.
            # qCovarColList is omitted when every categorical covariate was dropped
            # (sex is constant on a single-sex phenotype).
            cmd=(
                step1_fitNULLGLMM.R
                --bedFile="~{plink_bed_file}"
                --bimFile="~{plink_bim_file}"
                --famFile="~{plink_fam_file}"
                --useSparseGRMtoFitNULL=TRUE
                --useSparseGRMforVarRatio=TRUE
                --sparseGRMFile="~{sparseGRM}"
                --sparseGRMSampleIDFile="~{sparseGRM_IDlist}"
                --phenoFile="${pheno_file}"
                --phenoCol="${phecode}"
                --covarColList="${covar_for[${phecode}]}"
            )
            if [[ -n "${qcovar_for[${phecode}]}" ]]; then
                cmd+=(--qCovarColList="${qcovar_for[${phecode}]}")
            fi
            cmd+=(
                --sampleIDColinphenoFile=person_id
                --traitType=~{trait_type}
                --isCateVarianceRatio=TRUE
                --cateVarRatioMinMACVecExclude=~{cate_var_ratio_min_mac_exclude}
                --cateVarRatioMaxMACVecInclude=~{cate_var_ratio_max_mac_include}
                --IsOverwriteVarianceRatioFile=TRUE
                --outputPrefix="~{output_prefix}_step1Out_${phecode}"
            )
            if [[ "${have_include_opt}" -eq 1 && -n "${include_for[${phecode}]:-}" ]]; then
                cmd+=(--SampleIDIncludeFile="${include_for[${phecode}]}")
            fi
            if [[ -n "~{inv_norm_arg}" ]]; then
                cmd+=(~{inv_norm_arg})
            fi
            "${cmd[@]}"
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
        # phecode / status / covariates kept / categorical kept / dropped /
        # include file / per-phecode phenotype file. The last two are non-empty
        # exactly for the phecodes fit on one sex — check this before trusting
        # a lambda for a sex-specific phecode.
        File        covariate_plan          = "covariates_by_phecode.tsv"
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
        # false → omit --varianceRatioFile even if variance_ratios is set.
        Boolean use_varianceRatio = true
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
    # Same rule as SAIGE_phewas_Step2.wdl, additionally gated by use_varianceRatio.
    Boolean is_fastTest = use_varianceRatio && defined(variance_ratios)
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
        # Optional VR: use_varianceRatio=false, or undefined/empty files → omit
        # --varianceRatioFile (ratio=1). ${variance_ratio:+--varianceRatioFile=...}
        # below drops the flag when variance_ratio is empty.
        variance_ratios=(~{sep=' ' select_first([variance_ratios, []])})
        model_files=(~{sep=' ' GMMATmodelFiles})
        chromo_list=(~{sep=' ' chromo_list})
        use_vr=~{true="1" false="0" use_varianceRatio}
        have_vr=0
        if [[ "${use_vr}" -eq 1 && ${#variance_ratios[@]} -gt 0 ]]; then
            have_vr=1
        else
            echo "WARNING: Step 2 will omit --varianceRatioFile (use_varianceRatio=~{use_varianceRatio}, is_fastTest=~{is_fastTest})" >&2
        fi

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
        if [[ ${#model_files[@]} -eq 0 ]]; then
            echo "ERROR: Step 1 produced no null-model .rda files" >&2
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
            if [[ "${have_vr}" -eq 1 ]]; then
                for f in "${variance_ratios[@]}"; do
                    if [[ "$(basename "${f}")" == "${expected_vr}" && -s "${f}" ]]; then
                        variance_ratio="${f}"
                        break
                    fi
                done
            fi
            if [[ -z "${model_file}" ]]; then
                echo "WARNING: no null model matching ${phecode}; skipping" >&2
                echo "  expected rda: ${expected_rda}" >&2
                printf '%s\t%s\t%s\n' "${phecode}" "-" "missing_step1" >> failed_step2_singlevar.txt
                n_fail=$((n_fail + 1))
                continue
            fi
            if [[ "${have_vr}" -eq 1 && -z "${variance_ratio}" ]]; then
                echo "WARNING: no variance-ratio file matching ${phecode}; skipping" >&2
                echo "  expected vr:  ${expected_vr}" >&2
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
        # false → omit --varianceRatioFile even if variance_ratios is set.
        Boolean use_varianceRatio = true
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
        # Optional VR: use_varianceRatio=false, or undefined/empty files → omit
        # --varianceRatioFile (ratio=1).
        variance_ratios=(~{sep=' ' select_first([variance_ratios, []])})
        model_files=(~{sep=' ' GMMATmodelFiles})
        chromo_list=(~{sep=' ' chromo_list})
        use_vr=~{true="1" false="0" use_varianceRatio}
        have_vr=0
        if [[ "${use_vr}" -eq 1 && ${#variance_ratios[@]} -gt 0 ]]; then
            have_vr=1
        else
            echo "WARNING: Step 2 geneset will omit --varianceRatioFile (use_varianceRatio=~{use_varianceRatio})" >&2
        fi

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
        if [[ ${#model_files[@]} -eq 0 ]]; then
            echo "ERROR: Step 1 produced no null-model .rda files" >&2
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
            if [[ "${have_vr}" -eq 1 ]]; then
                for f in "${variance_ratios[@]}"; do
                    if [[ "$(basename "${f}")" == "${expected_vr}" && -s "${f}" ]]; then
                        variance_ratio="${f}"
                        break
                    fi
                done
            fi
            if [[ -z "${model_file}" ]]; then
                echo "WARNING: no null model matching ${phecode}; skipping" >&2
                echo "  expected rda: ${expected_rda}" >&2
                printf '%s\t%s\t%s\n' "${phecode}" "-" "missing_step1" >> failed_step2_geneset.txt
                n_fail=$((n_fail + 1))
                continue
            fi
            if [[ "${have_vr}" -eq 1 && -z "${variance_ratio}" ]]; then
                echo "WARNING: no variance-ratio file matching ${phecode}; skipping" >&2
                echo "  expected vr:  ${expected_vr}" >&2
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
