#!/usr/bin/env bash
set -euo pipefail

if (( $# < 2 )); then
    echo "Usage: seqkit_stats.sh GROUP_LABEL INPUT.fastq.gz [INPUT.fastq.gz ...]" >&2
    exit 2
fi

label=$1
shift
for input in "$@"; do
    if [[ ! -f $input || ! -r $input || $input != *.gz ]]; then
        echo "Expected a readable gzip FASTQ: $input" >&2
        exit 2
    fi
done

# Stream the FASTQ records as one group so Nx and quality are computed from
# all reads together, including both Illumina mates or all ONT chunks.
# Redirect stdout to a new TSV file to save the header and one result row.
gzip -cd -- "$@" | seqkit stats -a -N 90 -T -j 2 --stdin-label "$label" -
