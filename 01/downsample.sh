#!/usr/bin/env bash
set -euo pipefail
set -o noclobber

if [[ $# -lt 1 ]]; then
    echo "Usage: downsample.sh illumina TARGET_BP SEED R1 R2 OUT_R1 OUT_R2" >&2
    echo "       downsample.sh {ont|hifi} TARGET_BP SEED INPUT OUT" >&2
    exit 2
fi

platform=$1
shift
case "$platform:$#" in
    illumina:6|ont:4|hifi:4) ;;
    *) echo "Invalid platform or argument count. See README.md." >&2; exit 2 ;;
esac

target=$1
seed=$2
shift 2
[[ $target =~ ^[1-9][0-9]*$ && $seed =~ ^[0-9]+$ ]] || {
    echo "TARGET_BP must be a positive integer and SEED a nonnegative integer." >&2
    exit 2
}

# Use one combined base count for a paired library and one count for a long-read file.
if [[ $platform == illumina ]]; then
    inputs=("$1" "$2")
    outputs=("$3" "$4")
else
    inputs=("$1")
    outputs=("$2")
fi

for input in "${inputs[@]}"; do
    [[ -r $input && -f $input ]] || { echo "Cannot read input: $input" >&2; exit 2; }
done

for output in "${outputs[@]}"; do
    [[ ! -e $output && ! -L $output ]] || {
        echo "Output already exists: $output" >&2
        exit 2
    }
done
if [[ $platform == illumina && ${outputs[0]} == "${outputs[1]}" ]]; then
    echo "R1 and R2 outputs must be distinct." >&2
    exit 2
fi

# SeqKit's machine-readable summary reports total bases in the sum_len column.
input_bases=$(seqkit stats -T "${inputs[@]}" | awk -F '\t' '
    NR == 1 { for (i = 1; i <= NF; i++) if ($i == "sum_len") col = i; next }
    col && $col ~ /^[0-9]+$/ { total += $col }
    END { if (total <= 0) exit 1; printf "%.0f\n", total }
')
fraction=$(LC_ALL=C awk -v target="$target" -v total="$input_bases" '
    BEGIN { printf "%.17g\n", target >= total ? 1 : target / total }
')

# Seqtk samples whole reads. Identical seeds preserve mates only when the
# input R1 and R2 already contain the same pairs in the same order.
for index in "${!inputs[@]}"; do
    if [[ $fraction == 1 ]]; then
        seqtk seq "${inputs[index]}" | gzip -c > "${outputs[index]}"
    else
        seqtk sample -s"$seed" "${inputs[index]}" "$fraction" |
            gzip -c > "${outputs[index]}"
    fi
done

# Record actual output bases and read counts for comparison with the target.
seqkit stats -T "${outputs[@]}"
