#!/usr/bin/env bash
# Build a PLINK bed dataset from GRM VCFs and a SAIGE sparse GRM.
#
# LOW-DISK VARIANT of build_saige_plink_and_grm.sh. The GRM it produces is
# IDENTICAL -- same inputs, same filters, same createSparseGRM call. The only
# differences are that intermediates are deleted once consumed and the
# reheader copy is skipped when it is not needed. Both are disk optimisations
# and neither touches a single genotype.
#
# Why it exists: the original holds the localized VCFs, the bcftools concat
# output AND the bcftools reheader copy on disk simultaneously, so peak disk is
# about 3x the total input size before PLINK even starts. On 22 whole-chromosome
# FLARE anc VCFs (~100 GB total) that is ~300 GB of VCF alone, and the
# 2026-10-10 genome-wide run died with "No space left on device" during
# plink2 --vcf on a 300 GB disk. Freeing as we go takes peak to ~200 GB.
#
# Kept as a separate file rather than an edit: the build script is a File input
# to BuildSaigePlinkAndSparseGRM, so its hash is part of the call-caching key.
# Editing the original in place would invalidate the cached chr1+chr22 GRM and
# every null fit downstream of it.
#
# Inputs (env or flags):
#   --analysis-samples  sample ID list (VCF / analysis order)
#   --grm-vcfs          space-separated VCF paths (or pass as trailing args after --)
#   --relatedness-cutoff  SAIGE relatednessCutoff (default 0.05)
#   --n-threads
#   --out-prefix        output prefix directory/prefix (default: saige_grm/plink)
#
# Writes:
#   ${prefix}.bed/.bim/.fam
#   sparse GRM mtx + sampleIDs from createSparseGRM.R
#   saige_sparse_grm.paths.tsv documenting the produced file names

set -euo pipefail

ANALYSIS_SAMPLES=""
RELATEDNESS_CUTOFF="0.05"
N_THREADS="8"
OUT_PREFIX="saige_grm/plink"
NUM_RANDOM_MARKERS="2000"
CREATE_SPARSE_GRM_R="${CREATE_SPARSE_GRM_R:-/opt/SAIGE/extdata/createSparseGRM.R}"
VCFS=()

usage() {
  cat <<EOF
Usage: $0 --analysis-samples FILE --out-prefix PREFIX [options] -- vcf1 [vcf2 ...]
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --analysis-samples) ANALYSIS_SAMPLES="$2"; shift 2 ;;
    --relatedness-cutoff) RELATEDNESS_CUTOFF="$2"; shift 2 ;;
    --n-threads) N_THREADS="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    --num-random-markers) NUM_RANDOM_MARKERS="$2"; shift 2 ;;
    --create-sparse-grm-r) CREATE_SPARSE_GRM_R="$2"; shift 2 ;;
    --) shift; VCFS+=("$@"); break ;;
    -h|--help) usage; exit 0 ;;
    *)
      if [[ -f "$1" ]]; then
        VCFS+=("$1"); shift
      else
        echo "Unknown arg: $1" >&2; usage >&2; exit 1
      fi
      ;;
  esac
done

if [[ -z "${ANALYSIS_SAMPLES}" || ${#VCFS[@]} -eq 0 ]]; then
  echo "Required: --analysis-samples and at least one VCF" >&2
  usage >&2
  exit 1
fi
if [[ ! -x "${CREATE_SPARSE_GRM_R}" && ! -f "${CREATE_SPARSE_GRM_R}" ]]; then
  echo "createSparseGRM.R not found: ${CREATE_SPARSE_GRM_R}" >&2
  exit 1
fi

say_disk() {
  echo "[disk] $1: $(df -Pk . | awk 'NR==2{printf "%.1f GiB used, %.1f GiB free", $3/1048576, $4/1048576}')"
}

OUT_DIR="$(dirname "${OUT_PREFIX}")"
mkdir -p "${OUT_DIR}" grm_build
say_disk "start"
VCF_LIST=grm_build/vcf_list.txt
: > "${VCF_LIST}"
for vcf in "${VCFS[@]}"; do
  echo "${vcf}" >> "${VCF_LIST}"
done

N_VCF=$(wc -l < "${VCF_LIST}" | tr -d ' ')
if [[ "${N_VCF}" -eq 1 ]]; then
  INPUT_VCF=$(head -n1 "${VCF_LIST}")
else
  bcftools concat -f "${VCF_LIST}" -Oz -o grm_build/concat.vcf.gz --threads "${N_THREADS}"
  bcftools index -t grm_build/concat.vcf.gz
  INPUT_VCF=grm_build/concat.vcf.gz
  say_disk "after concat"
  # The concat is self-contained, so the localized inputs are dead weight.
  # They live on this task's own disk; Cromwell does not re-read them.
  for v in "${VCFS[@]}"; do rm -f "${v}" "${v}.tbi" "${v}.csi"; done
  CONCAT_IS_OURS=1
  say_disk "after freeing localized inputs"
fi

# FLARE anc VCFs often declare ##FORMAT=<ID=GT,...> twice; PLINK2 rejects that.
# Rewriting the file costs a full extra copy, so only do it when there IS a
# duplicate. Behaviour is unchanged either way: with one GT line the awk filter
# is a no-op and the reheader would reproduce the input byte for byte.
N_GT=$(bcftools view -h "${INPUT_VCF}" | grep -c '^##FORMAT=<ID=GT,' || true)
echo "FORMAT/GT header lines: ${N_GT}"
if [[ "${N_GT}" -gt 1 ]]; then
  bcftools view -h "${INPUT_VCF}" \
    | awk '/^##FORMAT=<ID=GT,/{ if (gt++) next } { print }' \
    > grm_build/header.fixed.txt
  bcftools reheader -h grm_build/header.fixed.txt -o grm_build/plink_in.vcf.gz "${INPUT_VCF}"
  bcftools index -t grm_build/plink_in.vcf.gz || bcftools index -c grm_build/plink_in.vcf.gz
  if [[ "${CONCAT_IS_OURS:-0}" -eq 1 ]]; then
    rm -f grm_build/concat.vcf.gz grm_build/concat.vcf.gz.tbi
    CONCAT_IS_OURS=0
  fi
  INPUT_VCF=grm_build/plink_in.vcf.gz
  PLINK_IN_IS_OURS=1
  say_disk "after reheader"
else
  echo "Single GT header line: skipping the reheader copy."
fi

plink2 --vcf "${INPUT_VCF}" \
  --double-id \
  --make-bed \
  --set-all-var-ids "@:#:\$r:\$a" \
  --out grm_build/all \
  --threads "${N_THREADS}"
say_disk "after plink2 import"
if [[ "${PLINK_IN_IS_OURS:-0}" -eq 1 ]]; then
  rm -f grm_build/plink_in.vcf.gz grm_build/plink_in.vcf.gz.tbi grm_build/plink_in.vcf.gz.csi
fi
if [[ "${CONCAT_IS_OURS:-0}" -eq 1 ]]; then
  rm -f grm_build/concat.vcf.gz grm_build/concat.vcf.gz.tbi
fi
say_disk "after freeing VCFs"

KEEP_PY="${MAKE_PLINK_KEEP_PY:-$(dirname "${BASH_SOURCE[0]}")/make_plink_keep.py}"
if [[ ! -f "${KEEP_PY}" ]]; then
  KEEP_PY="/opt/tractor_mix_scripts/make_plink_keep.py"
fi
python3 "${KEEP_PY}" \
  --analysis-samples "${ANALYSIS_SAMPLES}" \
  --fam grm_build/all.fam \
  --out-keep grm_build/keep_iid.txt \
  --out-samples grm_build/analysis_samples.intersect.txt

plink2 --bfile grm_build/all \
  --keep grm_build/keep_iid.txt \
  --make-bed \
  --out "${OUT_PREFIX}" \
  --threads "${N_THREADS}"
rm -f grm_build/all.bed grm_build/all.bim grm_build/all.fam
say_disk "before createSparseGRM"
echo "Markers entering createSparseGRM: $(wc -l < "${OUT_PREFIX}.bim" | tr -d ' ')"

# SAIGE createSparseGRM expects a PLINK prefix without extension
Rscript "${CREATE_SPARSE_GRM_R}" \
  --plinkFile="${OUT_PREFIX}" \
  --nThreads="${N_THREADS}" \
  --outputPrefix="${OUT_DIR}/sparseGRM" \
  --numRandomMarkerforSparseKin="${NUM_RANDOM_MARKERS}" \
  --relatednessCutoff="${RELATEDNESS_CUTOFF}"

# Locate outputs (SAIGE embeds cutoff + marker count in the filename)
SPARSE_MTX=$(ls -1 "${OUT_DIR}"/sparseGRM*.sparseGRM.mtx | head -n1)
SPARSE_IDS=$(ls -1 "${OUT_DIR}"/sparseGRM*.sparseGRM.mtx.sampleIDs.txt | head -n1)
if [[ -z "${SPARSE_MTX}" || -z "${SPARSE_IDS}" ]]; then
  echo "Failed to locate createSparseGRM outputs under ${OUT_DIR}" >&2
  ls -la "${OUT_DIR}" >&2 || true
  exit 1
fi

# Stable symlinks for WDL outputs
ln -sf "$(basename "${SPARSE_MTX}")" "${OUT_DIR}/sparseGRM.mtx"
ln -sf "$(basename "${SPARSE_IDS}")" "${OUT_DIR}/sparseGRM.sampleIDs.txt"

cat > "${OUT_DIR}/saige_sparse_grm.paths.tsv" <<EOF
key	path
plink_prefix	${OUT_PREFIX}
sparse_grm_mtx	${SPARSE_MTX}
sparse_grm_sample_ids	${SPARSE_IDS}
relatedness_cutoff	${RELATEDNESS_CUTOFF}
num_random_markers	${NUM_RANDOM_MARKERS}
EOF

echo "Wrote PLINK prefix ${OUT_PREFIX} and sparse GRM:"
echo "  ${SPARSE_MTX}"
echo "  ${SPARSE_IDS}"