#!/usr/bin/env bash
set -euo pipefail

if [[ $# != 5 || ! $5 =~ ^[1-9][0-9]*$ ]]; then
    echo "Usage: polypolish.sh DRAFT.fa[.gz] R1.fq.gz R2.fq.gz NEW_OUTPUT_DIR BWA_THREADS" >&2
    exit 2
fi

source_draft=$1
r1=$2
r2=$3
out=$4
threads=$5

for input in "$source_draft" "$r1" "$r2"; do
    if [[ ! -f $input || ! -r $input ]]; then
        echo "Cannot read input file: $input" >&2
        exit 2
    fi
done

# Create a new directory and index a local draft copy.
mkdir -- "$out"
draft="$out/draft.fasta"
if [[ $source_draft == *.gz ]]; then
    gzip -cd -- "$source_draft" > "$draft"
else
    cp -- "$source_draft" "$draft"
fi
bwa index "$draft" 2> "$out/bwa_index.log"

# Align each mate separately and retain all reported alignments.
bwa mem -a -t "$threads" "$draft" "$r1" > "$out/alignments_1.sam" 2> "$out/bwa_R1.log"
bwa mem -a -t "$threads" "$draft" "$r2" > "$out/alignments_2.sam" 2> "$out/bwa_R2.log"

# Use paired insert-size evidence to filter the two alignment files.
polypolish filter --in1 "$out/alignments_1.sam" --in2 "$out/alignments_2.sam" \
    --out1 "$out/filtered_1.sam" --out2 "$out/filtered_2.sam" \
    > "$out/polypolish_filter.log" 2>&1

# Run one conservative polishing round and verify FASTA identifiers and order.
polypolish polish --careful "$draft" "$out/filtered_1.sam" "$out/filtered_2.sam" \
    > "$out/final_contigs.fasta" 2> "$out/polypolish_polish.log"
test -s "$out/final_contigs.fasta"
diff -u <(awk '/^>/ {print $1}' "$draft") <(awk '/^>/ {print $1}' "$out/final_contigs.fasta")
gzip -c "$out/final_contigs.fasta" > "$out/contigs.fasta.gz"
