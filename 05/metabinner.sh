#!/usr/bin/env bash
set -euo pipefail

if [[ $# != 6 ]]; then
    echo "Usage: metabinner.sh METABINNER_DIR CONTIGS.fa[.gz] R1 R2 NEW_OUTPUT_DIR GROUP_ID" >&2
    exit 2
fi
root=$1
source=$2
r1=$3
r2=$4
out=$5
group=$6
helper="$(cd "$(dirname "$0")" && pwd)/profiles_and_bins.py"

# The official coverage script matches mates by their filename suffixes.
left=${r1##*/}
right=${r2##*/}
if [[ $r1 != *_1.cleaned.fastq.gz || $r2 != *_2.cleaned.fastq.gz ||
      ${left%_1.cleaned.fastq.gz} != "${right%_2.cleaned.fastq.gz}" ]]; then
    echo "Expected matching _1.cleaned.fastq.gz and _2.cleaned.fastq.gz files." >&2
    exit 2
fi

# Stage the reference under the basename used by the Workbench job.
mkdir -- "$out"
draft="$out/source_contigs.fasta"
filtered="$out/contigs_filtered.fasta"
if [[ $source == *.gz ]]; then
    gzip -cd -- "$source" > "$draft"
else
    cp -- "$source" "$draft"
fi

python3 "$helper" check-fasta "$draft"
python3 "$root/scripts/Filter_tooshort.py" "$draft" 1000 > "$out/filter_tooshort.log" 2>&1
cp -- "$out/source_contigs_1000.fa" "$filtered"
python3 "$helper" check-filter "$draft" "$filtered"
python3 "$root/scripts/gen_kmer.py" "$filtered" 1000 4 > "$out/gen_kmer.log" 2>&1

# Generate coverage using the same sample's cleaned Illumina mates.
coverage="$out/coverage_output"
bash "$root/scripts/gen_coverage_file.sh" -a "$filtered" -o "$coverage" \
    -t 16 -m 64 -l 1000 -f _1.cleaned.fastq.gz -r _2.cleaned.fastq.gz \
    "$r1" "$r2" > "$out/gen_coverage.log" 2>&1

# Reject mismatched IDs and put profile rows in FASTA order.
python3 "$helper" align "$filtered" "$coverage/coverage_profile.tsv" \
    "$out/contigs_filtered_kmer_4_f1000.csv" "$out"

target="$out/metabinner_output"
bash "$root/run_metabinner.sh" -a "$filtered" -o "$target" \
    -d "$out/coverage_profile.tsv" -k "$out/kmer_profile.csv" \
    -p "$root" -s large -t 8 > "$out/metabinner.log" 2>&1

# Export one FASTA for each final bin, retaining the unbinned contig list.
assignments="$target/metabinner_res/metabinner_result.tsv"
python3 "$helper" export "$filtered" "$assignments" "$out" "$group"
cp -- "$assignments" "$out/metabinner_result.tsv"
tar -czf "$out/raw_MAGs.tar.gz" -C "$out" raw_MAGs
