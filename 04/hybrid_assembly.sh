#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: hybrid_assembly.sh {ont|hifi} R1.fastq.gz R2.fastq.gz LONG.fastq.gz NEW_OUTPUT_DIR" >&2
    echo "       hybrid_assembly.sh resume EXISTING_OUTPUT_DIR" >&2
    exit 2
}

case "${1:-}:$#" in
    ont:5|hifi:5)
        # Use the long-read option corresponding to the input library.
        # The HiFi mode is an example and was not verified in the project runs.
        long_flag=--nanopore
        if [[ $1 == hifi ]]; then long_flag=--pacbio; fi
        for input in "$2" "$3" "$4"; do
            [[ -f $input && -r $input ]] || { echo "Unreadable FASTQ: $input" >&2; exit 2; }
        done
        out=$5
        [[ ! -e $out && ! -L $out ]] || { echo "Output already exists: $out" >&2; exit 2; }
        spades.py --meta --only-assembler -1 "$2" -2 "$3" "$long_flag" "$4" \
            -t 32 -m 256 --checkpoints all -o "$out"
        ;;
    resume:2)
        out=$2
        # Preserve SPAdes parameters and restart from the latest native checkpoint.
        if [[ ! -s $out/params.txt || ! -d $out/pipeline_state ]]; then
            echo "Missing SPAdes parameters or checkpoint state: $out" >&2
            exit 2
        fi
        OMP_NUM_THREADS=32 OMP_THREAD_LIMIT=32 OMP_DYNAMIC=FALSE \
            spades.py --continue -o "$out"
        ;;
    *) usage ;;
esac

# Confirm that SPAdes produced a nonempty final contig FASTA.
test -s "$out/contigs.fasta"
