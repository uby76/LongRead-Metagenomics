#!/usr/bin/env bash
set -euo pipefail
set -o noclobber

if [[ $# != 3 ]]; then
    echo "Usage: assemble_illumina.sh R1.fastq.gz R2.fastq.gz NEW_OUTPUT_DIR" >&2
    exit 2
fi

r1=$1
r2=$2
out=$3

# Assemble the downsampled, deduplicated mate files with the original preset.
# The output directory must be new. The memory argument is expressed in bytes.
megahit -1 "$r1" -2 "$r2" --presets meta-large \
    --min-contig-len 300 -t 16 -m 412316860416 -o "$out"

# Export the compressed assembly while preserving the assembler's native files.
gzip -c "$out/final.contigs.fa" > "$out/contigs.fasta.gz"
