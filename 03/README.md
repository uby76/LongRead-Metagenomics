# 03. Polypolish of ONT and PacBio HiFi contigs

The Workbench workflow polished each existing long-read assembly with the clean Illumina paired-end reads from the same sample. It used BWA 0.7.19 and Polypolish 0.7.1. The recovery run used the GNU/glibc Bioconda build of Polypolish. Both assembly types used the same alignment, filtering and polishing settings.

## Inputs

| Strategy | Draft assembly | Short reads |
| --- | --- | --- |
| ONT-polished | Existing ONT-only contigs from the cleaned, downsampled ONT library | Clean, downsampled Illumina R1 and R2 from the same sample |
| HiFi-polished | Existing HiFi-only contigs from the cleaned, downsampled HiFi library | The same clean, downsampled Illumina R1 and R2 |

Supply synchronized R1 and R2 files with matching read IDs in the same order. The Workbench worker checked the first 1,000 read pairs for matching IDs and valid FASTQ records before alignment.

## Commands used in Workbench

The worker first decompressed the draft to a local `draft.fasta` and indexed that copy. Assign `BWA_THREADS` to the number of CPUs allocated to the task. The commands below use `XXX` as the directory placeholder.

```bash
bwa index XXX/draft.fasta
bwa mem -a -t "$BWA_THREADS" XXX/draft.fasta XXX/sample_i_ds_1.cleaned.fastq.gz \
    > XXX/alignments_1.sam
bwa mem -a -t "$BWA_THREADS" XXX/draft.fasta XXX/sample_i_ds_2.cleaned.fastq.gz \
    > XXX/alignments_2.sam

polypolish filter --in1 XXX/alignments_1.sam --in2 XXX/alignments_2.sam \
    --out1 XXX/filtered_1.sam --out2 XXX/filtered_2.sam

polypolish polish --careful XXX/draft.fasta XXX/filtered_1.sam XXX/filtered_2.sam \
    > XXX/final_contigs.fasta
```

The separate BWA MEM calls use `-a` to retain multiple alignment placements. `polypolish filter` uses paired-end insert-size evidence to remove inconsistent placements. The workflow ran one `polypolish polish --careful` round using the default thresholds. Careful mode excludes reads with multiple remaining alignments from correction.

## Run the wrapper

Run from the repository root with Bash, BWA, Polypolish and gzip available. Provide a fresh output directory with an existing parent directory. The 16-thread allocation below is an example.

```bash
# ONT-polished
bash 03/polypolish.sh XXX/sample_n_ds.contigs.fasta.gz \
    XXX/sample_i_ds_1.cleaned.fastq.gz XXX/sample_i_ds_2.cleaned.fastq.gz \
    XXX/ont_polished 16

# HiFi-polished
bash 03/polypolish.sh XXX/sample_p_ds.contigs.fasta.gz \
    XXX/sample_i_ds_1.cleaned.fastq.gz XXX/sample_i_ds_2.cleaned.fastq.gz \
    XXX/hifi_polished 16
```

The wrapper retains the working draft, BWA index, original and filtered SAM files, and command logs in the output directory. It checks the final contig IDs and order against the draft and writes both `final_contigs.fasta` and `contigs.fasta.gz`. The historical worker also saved per-contig before and after lengths in `polishing_lengths.csv`.

[Polypolish running guide](https://github.com/rrwick/Polypolish/wiki/How-to-run-Polypolish)
