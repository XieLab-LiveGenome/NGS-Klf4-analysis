#!/usr/bin/env bash
set -euo pipefail

# Usage:
# ./mcool_to_fanc_ICE_VCSQRT_parallel.sh <input_mcool_folder> [output_root_folder]
#
# Example:
# nohup ./mcool_to_fanc_ICE_VCSQRT_parallel.sh ./cooler_output ./FANC_5kb_matrices_chr4 \
#     > fanc_5kb_ICE_VCSQRT_chr4.log 2>&1 &

INDIR="${1:-}"
OUTROOT="${2:-./FANC_5kb_matrices_chr4}"

if [[ -z "$INDIR" ]]; then
    echo "Usage: $0 <input_mcool_folder> [output_root_folder]"
    exit 1
fi

if [[ ! -d "$INDIR" ]]; then
    echo "Error: input folder does not exist: $INDIR"
    exit 1
fi

# -------------------------
# User settings
# -------------------------

FANC_ENV="/data2/xiongl2/venv_fanc/bin/activate"

RESOLUTION="5kb"
CHROM="chr4"

# Number of samples processed in parallel
JOBS=2

# Your original command used --downsample 0.99.
# If you truly want no downsampling, set DOWNSAMPLE=1
# or remove --downsample "$DOWNSAMPLE" from the fanc hic command.
DOWNSAMPLE=0.99

# Output folders
RAW_DIR="${OUTROOT}/RAW"
ICE_DIR="${OUTROOT}/ICE"
VCSQRT_DIR="${OUTROOT}/VCSQRT"

mkdir -p "$RAW_DIR" "$ICE_DIR" "$VCSQRT_DIR"

process_one() {
    mcool="$1"

    source "$FANC_ENV"

    base=$(basename "$mcool")
    sample=${base%.mcool}

    echo "========================================"
    echo "Processing sample: ${sample}"
    echo "Input mcool: ${mcool}"
    echo "Output root: ${OUTROOT}"
    echo "========================================"

    raw_hic="${RAW_DIR}/${sample}_${RESOLUTION}.raw.hic"

    # -------------------------
    # 1. Convert mcool@5kb to FAN-C .hic
    # -------------------------
    fanc hic \
        --deepcopy \
        --downsample "$DOWNSAMPLE" \
        -f "${mcool}@${RESOLUTION}" \
        "$raw_hic"

    # =========================================================
    # ICE normalization
    # =========================================================

    ice_hic="${ICE_DIR}/ICE.${sample}_${RESOLUTION}.hic"

    ice_oe_txt="${ICE_DIR}/${sample}_${RESOLUTION}.ICE.OE.${CHROM}.txt"
    ice_observed_txt="${ICE_DIR}/${sample}_${RESOLUTION}.ICE.observed.${CHROM}.txt"

    echo "Running ICE normalization for ${sample}"

    fanc hic \
        -n \
        -m ICE \
        "$raw_hic" \
        "$ice_hic"

    fanc dump \
        -s "$CHROM" \
        -e \
        "$ice_hic" \
        "$ice_oe_txt"

    fanc dump \
        -s "$CHROM" \
        "$ice_hic" \
        "$ice_observed_txt"

    # =========================================================
    # VC-SQRT normalization
    # =========================================================

    vcsqrt_hic="${VCSQRT_DIR}/VCSQRT.${sample}_${RESOLUTION}.hic"

    vcsqrt_oe_txt="${VCSQRT_DIR}/${sample}_${RESOLUTION}.VCSQRT.OE.${CHROM}.txt"
    vcsqrt_observed_txt="${VCSQRT_DIR}/${sample}_${RESOLUTION}.VCSQRT.observed.${CHROM}.txt"

    echo "Running VC-SQRT normalization for ${sample}"

    fanc hic \
        -n \
        -m VC-SQRT \
        "$raw_hic" \
        "$vcsqrt_hic"

    fanc dump \
        -s "$CHROM" \
        -e \
        "$vcsqrt_hic" \
        "$vcsqrt_oe_txt"

    fanc dump \
        -s "$CHROM" \
        "$vcsqrt_hic" \
        "$vcsqrt_observed_txt"

    echo "Finished sample: ${sample}"
}

export -f process_one
export OUTROOT RAW_DIR ICE_DIR VCSQRT_DIR
export FANC_ENV RESOLUTION CHROM DOWNSAMPLE

find "$INDIR" -type f -name "*.mcool" -print0 \
    | xargs -0 -n 1 -P "$JOBS" bash -c 'process_one "$0"'