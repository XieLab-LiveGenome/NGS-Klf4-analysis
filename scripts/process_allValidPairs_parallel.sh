#!/usr/bin/env bash
set -euo pipefail

# -------------------------
# Usage
# -------------------------
# ./process_allValidPairs_parallel.sh /path/to/allValidPairs_folder
# ./process_allValidPairs_parallel.sh /path/to/allValidPairs_folder /path/to/output_folder

INDIR="${1:-}"
OUTDIR="${2:-./processed_pairs}"

if [[ -z "$INDIR" ]]; then
    echo "Usage: $0 <input_folder> [output_folder]"
    exit 1
fi

if [[ ! -d "$INDIR" ]]; then
    echo "Error: input folder does not exist: $INDIR"
    exit 1
fi

# -------------------------
# User settings
# -------------------------
CAP_CHR="chr4"
CAP_START=55200000
CAP_END=55700000

# Number of samples processed in parallel
JOBS=4

# Threads used by pairtools sort per sample
SORT_THREADS=4

mkdir -p "$OUTDIR"

process_one() {
    infile="$1"

    sample=$(basename "$infile")
    sample=${sample%.allValidPairs}

    echo "Processing ${sample}"

    pairs_gz="${OUTDIR}/${sample}.contacts.pairs.gz"
    header_gz="${OUTDIR}/${sample}.contacts.with_header.pairs.gz"
    sorted_gz="${OUTDIR}/${sample}.contacts.sorted.pairs.gz"
    selected_gz="${OUTDIR}/${sample}.capture_${CAP_CHR}_${CAP_START}_${CAP_END}.selected.pairs.gz"

    # 1. Convert HiC-Pro allValidPairs to pairtools pairs format
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5,$6,$4,$7}' "$infile" \
        | bgzip -c > "$pairs_gz"

    # 2. Add pairs header
    {
        echo "## pairs format v1.0"
        echo "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2"
        zcat "$pairs_gz"
    } | bgzip -c > "$header_gz"

    # 3. Sort pairs
    pairtools sort \
        --nproc "$SORT_THREADS" \
        -o "$sorted_gz" \
        "$header_gz"

    # 4. Select capture region
    pairtools select \
        "(chrom1 == '${CAP_CHR}' and pos1 >= ${CAP_START} and pos1 <= ${CAP_END}) and (chrom2 == '${CAP_CHR}' and pos2 >= ${CAP_START} and pos2 <= ${CAP_END})" \
        -o "$selected_gz" \
        "$sorted_gz"

    echo "Finished ${sample}"
}

export -f process_one
export OUTDIR CAP_CHR CAP_START CAP_END SORT_THREADS

find "$INDIR" -type f -name "*.allValidPairs" -print0 \
    | xargs -0 -n 1 -P "$JOBS" bash -c 'process_one "$0"'