# 01. Seqtk downsampling

For read counts, length and base quality metrics, see [SeqKit statistics](seqkit_stats.md).

This step downsamples Illumina paired-end, Oxford Nanopore (ONT), and PacBio HiFi FASTQ files to a sample-specific target number of bases. Each input and output path below uses `XXX` as a placeholder. Run the commands from the repository root.

## Project records

The project manifest assigned the same target to the Illumina, ONT and HiFi libraries of each sample. Across the seven samples, the targets were 23.1, 29.4, 30.8, 32.7, 33.6, 34.4 and 38.0 Gb. The ONT input manifest selected pass reads. The published downsampled files were two mate files for Illumina and one FASTQ for each long-read library.

The original Seqtk invocation and random seed could not be verified. The commands below are a reusable implementation. The example seed is 100, and the example target is 23.1 billion bases. Use the target recorded for the sample being processed.

## Requirements and input preparation

Install Bash, [Seqtk](https://github.com/lh3/seqtk), SeqKit and gzip. Prepare one complete FASTQ for each mate and one FASTQ per long-read library. Supply synchronized Illumina mates with the same read IDs in the same order. For ONT, supply the combined pass FASTQ for the sample. For PacBio, supply the HiFi FASTQ for the sample.

## Run one sample

```bash
# Illumina: the base target covers both mates combined.
bash 01/downsample.sh illumina 23100000000 100 \
    XXX/illumina_R1.fastq.gz XXX/illumina_R2.fastq.gz \
    XXX/illumina_R1.downsampled.fastq.gz XXX/illumina_R2.downsampled.fastq.gz

# ONT: sample whole pass reads.
bash 01/downsample.sh ont 23100000000 100 \
    XXX/ont_pass.fastq.gz XXX/ont.downsampled.fastq.gz

# PacBio HiFi: sample whole HiFi reads.
bash 01/downsample.sh hifi 23100000000 100 \
    XXX/hifi.fastq.gz XXX/hifi.downsampled.fastq.gz
```

The script divides the target bases by the input bases and passes that fraction to `seqtk sample`. It uses the same fraction and seed for both Illumina mates. When the input is already at or below the target, it retains all records. The reported base count can differ from the target because Seqtk samples reads rather than individual bases. The script prints observed read counts and base totals and refuses to overwrite output files.
