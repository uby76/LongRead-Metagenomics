#!/usr/bin/env bash
set -euo pipefail

if [[ $# != 3 ]]; then
    echo "Usage: assemble_long_reads.sh {ont|hifi} INPUT.fastq.gz NEW_OUTPUT_DIR" >&2
    exit 2
fi

case "$1" in
    ont) input_option=--in-ont ;;
    hifi) input_option=--in-hifi ;;
    *) echo "Platform must be ont or hifi." >&2; exit 2 ;;
esac

input=$2
out=$3
if [[ -e $out || -L $out ]]; then
    echo "Output path already exists: $out" >&2
    exit 2
fi

# Select the platform-specific MetaMDBG mode for one cleaned, downsampled library.
# The recorded commands used the default minimizer settings of MetaMDBG v1.4.
metaMDBG --out-dir "$out" "$input_option" "$input" \
    --threads 8 --min-contig-length 300
