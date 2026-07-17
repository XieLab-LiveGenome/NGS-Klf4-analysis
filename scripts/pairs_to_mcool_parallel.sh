#!/usr/bin/env bash
set -euo pipefail

# Usage:
# ./pairs_to_mcool_parallel.sh <input_folder> [output_folder]

INDIR="${1:-}"
OUTDIR="${2:-./cooler_output}"

if [[ -z "$INDIR" ]]; then
    echo "Usage: $0 <input_folder> [output_folder]"
    exit 1
fi

if [[ ! -d "$INDIR" ]]; then
    echo "Error: input folder does not exist: $INDIR"
    exit 1
fi

mkdir -p "$OUTDIR"

# -------------------------
# User settings
# -------------------------

CHROMSIZES="/data2/xiongl2/Softwares/HiC_PRO_INSTALL3.1.0/HiC-Pro_3.1.0/chrom_mm10.sizes"
ASSEMBLY="mm10"

# Initial resolution for cooler cload
BASE_RES=100

# Resolutions for mcool
RESOLUTIONS="200,500,1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000"

# Number of samples processed in parallel
JOBS=4

# Threads per sample for cooler zoomify
ZOOM_THREADS=10

# Input pattern
PATTERN="*.capture_chr4_55200000_55700000.selected.pairs.gz"

process_one() {
    infile="$1"

    base=$(basename "$infile")
    sample=${base%.capture_chr4_55200000_55700000.selected.pairs.gz}

    echo "Processing ${sample}"

    cool_file="${OUTDIR}/${sample}.cool"
    mcool_file="${OUTDIR}/${sample}.mcool"

    # 1. Create .cool from selected pairs
    cooler cload pairs \
        -c1 2 \
        -p1 3 \
        -c2 4 \
        -p2 5 \
        --assembly "$ASSEMBLY" \
        "${CHROMSIZES}:${BASE_RES}" \
        "$infile" \
        "$cool_file"

    # 2. Create balanced multi-resolution .mcool
    cooler zoomify \
        --nproc "$ZOOM_THREADS" \
        --out "$mcool_file" \
        --resolutions "$RESOLUTIONS" \
        --balance \
        "$cool_file"

    echo "Finished ${sample}"
}

export -f process_one
export OUTDIR CHROMSIZES ASSEMBLY BASE_RES RESOLUTIONS ZOOM_THREADS

find "$INDIR" -type f -name "$PATTERN" -print0 \
    | xargs -0 -n 1 -P "$JOBS" bash -c 'process_one "$0"'