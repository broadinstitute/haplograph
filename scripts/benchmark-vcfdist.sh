#!/usr/bin/env bash
# End-to-end haplograph benchmark: call → merge → vcfdist vs truth.
#
# Example:
#   ./scripts/benchmark-vcfdist.sh \
#     --bam         /path/to/HG002.bam \
#     --reference   /path/to/GRCh38.fa \
#     --locus       "chr6:28500000-33500000" \
#     --truth-vcf   /path/to/HG002.truth.vcf.gz \
#     --sample      HG002 \
#     --output-dir  ./benchmark/HG002_MHC

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

LOG_FILE=""
WORK_DIR=""

log() {
    # tee -a when logging is active; otherwise stderr only (early validation).
    if [[ -n "$LOG_FILE" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
    fi
}

# Re-compress plain-text VCFs that were written with a .vcf.gz suffix (legacy merge output).
ensure_bgzf_vcf() {
    local vcf="$1"
    if [[ ! -f "$vcf" ]]; then
        return 0
    fi
    if file -b "$vcf" | grep -qi 'gzip'; then
        return 0
    fi
    log "Re-compressing plain VCF to BGZF: ${vcf}"
    local tmp="${vcf}.bgzf.tmp"
    bcftools view -Oz -o "$tmp" -- "$vcf"
    mv -f "$tmp" "$vcf"
    rm -f "${vcf}.tbi" "${vcf}.csi"
}

# Copy a gs:// BAM to a local indexed BAM so parallel haplograph reads do not
# hammer the remote file (HTSlib GCS seeks fail under heavy concurrency).
# Sets LOCAL_BAM to the local path (stdout is not used — log() goes to stdout).
localize_remote_bam() {
    local remote_bam="$1"
    local locus="$2"
    local dest_dir="$3"
    local sample="$4"
    local locus_tag="$5"

    LOCAL_BAM="${dest_dir}/${sample}_${locus_tag}.local.bam"

    if [[ -f "$LOCAL_BAM" ]] && samtools quickcheck -v "$LOCAL_BAM" 2>/dev/null; then
        log "Using cached local BAM: ${LOCAL_BAM}"
    else
        if [[ -f "$LOCAL_BAM" ]]; then
            log "Cached local BAM failed quickcheck; re-downloading: ${LOCAL_BAM}"
            rm -f "$LOCAL_BAM" "${LOCAL_BAM}.bai"
        fi
        log "Localizing gs:// BAM for locus ${locus}"
        log "Command: samtools view -bhX ${remote_bam} ${remote_bam}.bai ${locus} > ${LOCAL_BAM}"
        samtools view -bhX "$remote_bam" "$locus" > "$LOCAL_BAM"
        log "Command: samtools index ${LOCAL_BAM}"
        samtools index "$LOCAL_BAM"
        log "Local BAM ready: ${LOCAL_BAM} ($(du -h "$LOCAL_BAM" | cut -f1))"
    fi

    if [[ ! -f "${LOCAL_BAM}.bai" ]]; then
        log "Command: samtools index ${LOCAL_BAM}"
        samtools index "$LOCAL_BAM"
    fi
}

on_error() {
    local exit_code=$?
    local line="${BASH_LINENO[0]:-unknown}"
    if [[ -n "$LOG_FILE" ]]; then
        {
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: failed at line ${line} (exit ${exit_code})"
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Log file: ${LOG_FILE}"
        } >>"$LOG_FILE"
    fi
    echo "ERROR: benchmark failed at line ${line} (exit ${exit_code}). See log: ${LOG_FILE:-stderr}" >&2
    exit "$exit_code"
}

on_exit() {
    local exit_code=$?
    if [[ $exit_code -ne 0 && -n "$LOG_FILE" ]]; then
        log "Benchmark exited with status ${exit_code}"
    fi
    if [[ -n "$WORK_DIR" && "$KEEP_TEMP" == false ]]; then
        rm -rf "$WORK_DIR"
    fi
}

trap on_error ERR
trap on_exit EXIT

usage() {
    cat <<'EOF'
Usage: benchmark-vcfdist.sh [options]

End-to-end pipeline:
  1. haplograph haplograph  (per sub-locus if locus > maximal-locus-size)
  2. haplograph dev-tools merge (skipped for single-chunk loci)
  3. subset truth VCF to the locus
  4. vcfdist query vs truth (merged germline VCF for full locus)

Required:
  --bam PATH              Input BAM (local path or gs://)
  --reference PATH        Reference FASTA
  --locus REGION          Genomic region chrom:start-end
  --truth-vcf PATH        Truth VCF
  --sample NAME           Sample ID (query VCF / haplograph --sampleid)

Optional:
  --truth-sample NAME     Truth VCF sample (default: same as --sample)
  --output-dir DIR        Run directory (default: ./benchmark_out)
  --haplograph-bin PATH   haplograph binary (default: $HAPLOGRAPH_BIN,
                          else target/release/haplograph, else PATH)
  --maximal-locus-size N  Sub-locus size for large regions (default: 200000)
  --overlap-bp N          Sub-locus overlap (default: 500)
  --min-mean-depth F      Skip sub-loci below this depth (default: 2.0)
  --window-size N         haplograph window size (default: 100)
  --coverage-tag STR      Label in vcfdist output filenames (default: NA)
  --verbosity N           vcfdist -v level (default: 1)
  --extra-vcfdist-args STR  Extra vcfdist arguments (quoted)
  --skip-haplograph       Skip steps 1–2; use existing outputs under output-dir
  --keep-temp             Keep intermediate VCFs
  --verbose               Pass --verbose to haplograph / merge
  -h, --help              Show this help

Requires: bcftools, samtools, vcfdist, and a built haplograph binary.

Remote gs:// BAMs are copied locally before haplograph runs:
  samtools view -bhX <gs://bam> <locus> > <output-dir>/haplograph/<sample>_<locus>.local.bam
(Reuses the cached file on subsequent runs if samtools quickcheck passes.)

Key outputs:
  <output-dir>/benchmark.log                          full run log (stdout + stderr)
  <output-dir>/haplograph/<prefix>.germline.vcf.gz   merged (or direct) callset
  <output-dir>/vcfdist/<prefix>.<locus>.<coverage>.precision-recall-summary.tsv
EOF
}

BAM=""
REFERENCE=""
LOCUS=""
TRUTH_VCF=""
SAMPLE=""
TRUTH_SAMPLE=""
OUTPUT_DIR="./benchmark_out"
HAPLOGRAPH_BIN="${HAPLOGRAPH_BIN:-}"
MAXIMAL_LOCUS_SIZE=200000
OVERLAP_BP=500
MIN_MEAN_DEPTH=2.0
WINDOW_SIZE=100
COVERAGE_TAG="NA"
VERBOSITY=1
EXTRA_VCFDIST_ARGS=""
SKIP_HAPLOGRAPH=false
KEEP_TEMP=false
VERBOSE=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam) BAM="$2"; shift 2 ;;
        --reference) REFERENCE="$2"; shift 2 ;;
        --locus) LOCUS="${2//,/}"; shift 2 ;;
        --truth-vcf) TRUTH_VCF="$2"; shift 2 ;;
        --sample) SAMPLE="$2"; shift 2 ;;
        --truth-sample) TRUTH_SAMPLE="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --haplograph-bin) HAPLOGRAPH_BIN="$2"; shift 2 ;;
        --maximal-locus-size) MAXIMAL_LOCUS_SIZE="$2"; shift 2 ;;
        --overlap-bp) OVERLAP_BP="$2"; shift 2 ;;
        --min-mean-depth) MIN_MEAN_DEPTH="$2"; shift 2 ;;
        --window-size) WINDOW_SIZE="$2"; shift 2 ;;
        --coverage-tag) COVERAGE_TAG="$2"; shift 2 ;;
        --verbosity) VERBOSITY="$2"; shift 2 ;;
        --extra-vcfdist-args) EXTRA_VCFDIST_ARGS="$2"; shift 2 ;;
        --skip-haplograph) SKIP_HAPLOGRAPH=true; shift ;;
        --keep-temp) KEEP_TEMP=true; shift ;;
        --verbose) VERBOSE=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

resolve_vcfdist_bin() {
    if command -v vcfdist >/dev/null 2>&1; then
        command -v vcfdist
        return 0
    fi
    local conda_vcfdist="/opt/miniconda3/envs/vcfdist/bin/vcfdist"
    if [[ -x "$conda_vcfdist" ]]; then
        echo "$conda_vcfdist"
        return 0
    fi
    return 1
}

resolve_haplograph_bin() {
    if [[ -n "$HAPLOGRAPH_BIN" ]]; then
        echo "$HAPLOGRAPH_BIN"
        return
    fi
    if [[ -x "${REPO_ROOT}/target/release/haplograph" ]]; then
        echo "${REPO_ROOT}/target/release/haplograph"
        return
    fi
    if [[ -x "${REPO_ROOT}/target/debug/haplograph" ]]; then
        echo "${REPO_ROOT}/target/debug/haplograph"
        return
    fi
    command -v haplograph
}

missing=()
[[ -z "$BAM" ]] && missing+=("--bam")
[[ -z "$REFERENCE" ]] && missing+=("--reference")
[[ -z "$LOCUS" ]] && missing+=("--locus")
[[ -z "$TRUTH_VCF" ]] && missing+=("--truth-vcf")
[[ -z "$SAMPLE" ]] && missing+=("--sample")
if ((${#missing[@]} > 0)); then
    echo "Error: missing required option(s): ${missing[*]}" >&2
    usage >&2
    exit 1
fi

if [[ -z "$TRUTH_SAMPLE" ]]; then
    TRUTH_SAMPLE="$SAMPLE"
fi

if ! command -v bcftools >/dev/null 2>&1; then
    echo "Error: 'bcftools' not found on PATH." >&2
    exit 1
fi

if [[ "$BAM" == gs://* && "$SKIP_HAPLOGRAPH" == false ]] && ! command -v samtools >/dev/null 2>&1; then
    echo "Error: 'samtools' not found on PATH (required to localize gs:// BAMs)." >&2
    exit 1
fi

VCFDIST_BIN="$(resolve_vcfdist_bin)" || {
    echo "Error: 'vcfdist' not found on PATH (try: conda activate vcfdist)." >&2
    exit 1
}

HAPLOGRAPH_BIN="$(resolve_haplograph_bin)" || {
    echo "Error: haplograph binary not found. Build with 'cargo build --release' or set --haplograph-bin." >&2
    exit 1
}

if [[ "$SKIP_HAPLOGRAPH" == false ]]; then
    if [[ ! -f "$REFERENCE" ]]; then
        echo "Error: reference not found: $REFERENCE" >&2
        exit 1
    fi
    if [[ "$BAM" != gs://* && ! -f "$BAM" ]]; then
        echo "Error: BAM not found: $BAM" >&2
        exit 1
    fi
fi

if [[ ! -f "$TRUTH_VCF" ]]; then
    echo "Error: truth VCF not found: $TRUTH_VCF" >&2
    exit 1
fi

REFERENCE="$(cd "$(dirname "$REFERENCE")" && pwd)/$(basename "$REFERENCE")"
TRUTH_VCF="$(cd "$(dirname "$TRUTH_VCF")" && pwd)/$(basename "$TRUTH_VCF")"

mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

LOG_FILE="${OUTPUT_DIR}/benchmark.log"
: >"$LOG_FILE"
exec > >(tee -a "$LOG_FILE") 2>&1

LOCUS_TAG="${LOCUS//:/_}"
LOCUS_TAG="${LOCUS_TAG//-/_}"
RUN_PREFIX="${OUTPUT_DIR}/haplograph/${SAMPLE}_${LOCUS_TAG}"
VCFDIST_DIR="${OUTPUT_DIR}/vcfdist"
HAPLOGRAPH_DIR="${OUTPUT_DIR}/haplograph"
mkdir -p "$HAPLOGRAPH_DIR" "$VCFDIST_DIR"

WORK_DIR="$(mktemp -d "${OUTPUT_DIR}/.vcfdist_work.XXXXXX")"
export TMPDIR="$WORK_DIR"

LOCUS_SPAN="${LOCUS#*:}"
LOCUS_START="${LOCUS_SPAN%-*}"
LOCUS_END="${LOCUS_SPAN#*-}"
LOCUS_LEN=$((LOCUS_END - LOCUS_START))
if (( LOCUS_LEN <= 0 )); then
    log "ERROR: invalid locus span (${LOCUS_START}-${LOCUS_END} => ${LOCUS_LEN} bp). Use digits only, e.g. chr6:29962395-30315226"
    exit 1
fi

CHUNKS=1
if (( LOCUS_LEN > MAXIMAL_LOCUS_SIZE )); then
    step=$((MAXIMAL_LOCUS_SIZE - OVERLAP_BP))
    CHUNKS=$(( (LOCUS_LEN + step - 1) / step ))
fi

log "Benchmark started"
log "  output-dir:          ${OUTPUT_DIR}"
log "  log-file:            ${LOG_FILE}"
log "  bam:                 ${BAM}"
log "  reference:           ${REFERENCE}"
log "  locus:               ${LOCUS} (${LOCUS_LEN} bp, ~${CHUNKS} sub-loci)"
log "  sample:              ${SAMPLE} (truth: ${TRUTH_SAMPLE})"
log "  haplograph:          ${HAPLOGRAPH_BIN}"
log "  vcfdist:             ${VCFDIST_BIN}"
log "  run-prefix:          ${RUN_PREFIX}"
log "  maximal-locus-size:  ${MAXIMAL_LOCUS_SIZE}"
log "  overlap-bp:          ${OVERLAP_BP}"
log "  min-mean-depth:      ${MIN_MEAN_DEPTH}"
log "  skip-haplograph:     ${SKIP_HAPLOGRAPH}"
if [[ "$SKIP_HAPLOGRAPH" == false && $CHUNKS -gt 1 ]]; then
    log "NOTE: large loci run haplograph per sub-locus; the first VCF appears only after sub-locus 0 finishes."
fi

HAPLOGRAPH_BAM="$BAM"
if [[ "$SKIP_HAPLOGRAPH" == false && "$BAM" == gs://* ]]; then
    log "[0/4] localize gs:// BAM"
    localize_remote_bam "$BAM" "$LOCUS" "$HAPLOGRAPH_DIR" "$SAMPLE" "$LOCUS_TAG"
    HAPLOGRAPH_BAM="$LOCAL_BAM"
    log "  localized BAM: ${HAPLOGRAPH_BAM}"
fi

if [[ "$SKIP_HAPLOGRAPH" == false ]]; then
    log "[1/4] haplograph haplograph"
    log "  BAM:      ${HAPLOGRAPH_BAM}"
    log "  Locus:    ${LOCUS} (${LOCUS_LEN} bp)"
    log "  Prefix:   ${RUN_PREFIX}"
    HAPLOGRAPH_ARGS=(
        --alignment-bam "$HAPLOGRAPH_BAM"
        --reference-fa "$REFERENCE"
        --sampleid "$SAMPLE"
        --output-prefix "$RUN_PREFIX"
        --locus "$LOCUS"
        --file-format vcf
        --window-size "$WINDOW_SIZE"
        --maximal-locus-size "$MAXIMAL_LOCUS_SIZE"
        --overlap-bp "$OVERLAP_BP"
        --min-mean-depth "$MIN_MEAN_DEPTH"
    )
    [[ "$VERBOSE" == true ]] && HAPLOGRAPH_ARGS+=(--verbose)
    log "Command: ${HAPLOGRAPH_BIN} haplograph ${HAPLOGRAPH_ARGS[*]}"
    "$HAPLOGRAPH_BIN" haplograph "${HAPLOGRAPH_ARGS[@]}"
    log "haplograph haplograph finished"

    if (( LOCUS_LEN > MAXIMAL_LOCUS_SIZE )); then
        log "[2/4] haplograph dev-tools merge (large locus)"
        MERGE_ARGS=(
            --output-prefix "$RUN_PREFIX"
            --locus "$LOCUS"
            --reference-fa "$REFERENCE"
            --sampleid "$SAMPLE"
            --maximal-locus-size "$MAXIMAL_LOCUS_SIZE"
            --overlap-bp "$OVERLAP_BP"
        )
        [[ "$VERBOSE" == true ]] && MERGE_ARGS+=(--verbose)
        log "Command: ${HAPLOGRAPH_BIN} dev-tools merge ${MERGE_ARGS[*]}"
        "$HAPLOGRAPH_BIN" dev-tools merge "${MERGE_ARGS[@]}"
        log "haplograph dev-tools merge finished"
    else
        log "[2/4] haplograph dev-tools merge skipped (locus fits in one sub-locus)"
    fi
else
    log "[1–2/4] haplograph steps skipped (--skip-haplograph)"
fi

QUERY_VCF="${RUN_PREFIX}.germline.vcf.gz"
if [[ ! -f "$QUERY_VCF" ]]; then
    log "ERROR: expected query VCF not found: ${QUERY_VCF}"
    if (( LOCUS_LEN > MAXIMAL_LOCUS_SIZE )); then
        log "Segment outputs present:"
        ls -la "${RUN_PREFIX}"_* 2>&1 || true
    else
        log "haplograph dir listing:"
        ls -la "$HAPLOGRAPH_DIR" 2>&1 || true
    fi
    exit 1
fi
QUERY_VCF="$(cd "$(dirname "$QUERY_VCF")" && pwd)/$(basename "$QUERY_VCF")"

log "[3/4] subset truth VCF to ${LOCUS} (sample ${TRUTH_SAMPLE})"
TRUTH_SUBSET="${WORK_DIR}/${SAMPLE}.${LOCUS_TAG}.base.vcf.gz"
if [[ ! -f "${TRUTH_VCF}.tbi" && ! -f "${TRUTH_VCF}.csi" ]]; then
    bcftools index -t "$TRUTH_VCF"
fi
bcftools view -s "$TRUTH_SAMPLE" -r "$LOCUS" "$TRUTH_VCF" -Oz -o "$TRUTH_SUBSET"
if [[ ! -f "$TRUTH_SUBSET" ]]; then
    log "ERROR: failed to create truth subset VCF: ${TRUTH_SUBSET}"
    exit 1
fi
bcftools index -t "$TRUTH_SUBSET"

QUERY_FILTERED="${WORK_DIR}/${SAMPLE}.filtered.query.vcf.gz"
log "[4/4] vcfdist (query: ${QUERY_VCF})"
ensure_bgzf_vcf "$QUERY_VCF"
if [[ ! -f "${QUERY_VCF}.tbi" && ! -f "${QUERY_VCF}.csi" ]]; then
    bcftools index -t "$QUERY_VCF"
fi
bcftools filter -e 'INFO/SOMATIC=1' "$QUERY_VCF" -Oz -o "$QUERY_FILTERED"
bcftools index -t "$QUERY_FILTERED"

cd "$WORK_DIR"
# shellcheck disable=SC2086
log "Command: ${VCFDIST_BIN} $(basename "$QUERY_FILTERED") $(basename "$TRUTH_SUBSET") ${REFERENCE} -v ${VERBOSITY} ${EXTRA_VCFDIST_ARGS}"
"$VCFDIST_BIN" \
    "$(basename "$QUERY_FILTERED")" \
    "$(basename "$TRUTH_SUBSET")" \
    "$REFERENCE" \
    -v "$VERBOSITY" \
    $EXTRA_VCFDIST_ARGS
log "vcfdist finished"

SUMMARY_VCF="summary.vcf"
if [[ ! -f "$SUMMARY_VCF" ]]; then
    log "ERROR: vcfdist did not produce summary.vcf in ${WORK_DIR}"
    ls -la "$WORK_DIR" 2>&1 || true
    exit 1
fi

OUT_TAG="${SAMPLE}.${LOCUS_TAG}.${COVERAGE_TAG}"
for tsv in *.tsv; do
    [[ -e "$tsv" ]] || continue
    mv "$tsv" "${VCFDIST_DIR}/${OUT_TAG}.${tsv}"
done
mv "$SUMMARY_VCF" "${VCFDIST_DIR}/${OUT_TAG}.summary.vcf"
cd - >/dev/null

PR_SUMMARY="${VCFDIST_DIR}/${OUT_TAG}.precision-recall-summary.tsv"
PHASING_SUMMARY="${VCFDIST_DIR}/${OUT_TAG}.phasing-summary.tsv"

log "Benchmark complete"
log "  query VCF:    ${QUERY_VCF}"
log "  truth subset: ${TRUTH_SUBSET}"
log "  vcfdist dir:  ${VCFDIST_DIR}"
log "  log file:     ${LOG_FILE}"
echo

if [[ -f "$PR_SUMMARY" ]]; then
    echo "==> precision-recall-summary.tsv"
    column -t -s $'\t' "$PR_SUMMARY" 2>/dev/null || cat "$PR_SUMMARY"
    echo
fi

if [[ -f "$PHASING_SUMMARY" ]]; then
    echo "==> phasing-summary.tsv"
    column -t -s $'\t' "$PHASING_SUMMARY" 2>/dev/null || cat "$PHASING_SUMMARY"
    echo
fi

if [[ "$KEEP_TEMP" == true ]]; then
    mkdir -p "${VCFDIST_DIR}/intermediate"
    cp -f "$TRUTH_SUBSET" "${VCFDIST_DIR}/intermediate/"
    cp -f "${TRUTH_SUBSET}.tbi" "${VCFDIST_DIR}/intermediate/" 2>/dev/null || true
    cp -f "$QUERY_FILTERED" "${VCFDIST_DIR}/intermediate/"
    cp -f "${QUERY_FILTERED}.tbi" "${VCFDIST_DIR}/intermediate/" 2>/dev/null || true
    echo "Intermediate VCFs: ${VCFDIST_DIR}/intermediate/"
fi
