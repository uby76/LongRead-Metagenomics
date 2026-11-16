# 01. Seqtk downsampling

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

## SeqKit statistics for raw and downsampled FASTQ

The original analysis calculated one row per dataset, sample and sequencing platform using `seqkit stats -a -N 90 -T -j 2`. It pooled each group through standard input and assigned it a label. Use Bash, gzip and a SeqKit release that supports `stats -N`. The examples below use gzip-compressed FASTQ and the default Sanger (Phred+33) quality encoding. Run them from the repository root. Replace each `XXX` with the appropriate directory.

### Raw data

```bash
# Illumina: include R1 and R2, plus any additional lanes for the same sample.
bash 01/seqkit_stats.sh raw_sample_illumina \
    XXX/sample_R1.fastq.gz XXX/sample_R2.fastq.gz > XXX/raw_sample_illumina.tsv

# ONT: include every pass FASTQ chunk for this sample.
bash 01/seqkit_stats.sh raw_sample_ont \
    XXX/ont_pass/*.fastq.gz > XXX/raw_sample_ont.tsv

# PacBio HiFi: use one FASTQ copy of the sample.
bash 01/seqkit_stats.sh raw_sample_hifi \
    XXX/sample_hifi.fastq.gz > XXX/raw_sample_hifi.tsv
```

The raw ONT scope is pass reads. The raw HiFi extended run used an existing FASTQ copy and cross-checked its read counts and length statistics against the earlier PacBio BAM index (PBI) results. Base-level quality and GC in this workflow come from FASTQ records.

### Downsampled data

```bash
# Illumina: pool the two downsampled mate files.
bash 01/seqkit_stats.sh downsampled_sample_illumina \
    XXX/sample_R1.downsampled.fastq.gz XXX/sample_R2.downsampled.fastq.gz \
    > XXX/downsampled_sample_illumina.tsv

# ONT: use the sample's downsampled ONT FASTQ.
bash 01/seqkit_stats.sh downsampled_sample_ont \
    XXX/sample_ont.downsampled.fastq.gz > XXX/downsampled_sample_ont.tsv

# PacBio HiFi: use the sample's downsampled HiFi FASTQ.
bash 01/seqkit_stats.sh downsampled_sample_hifi \
    XXX/sample_hifi.downsampled.fastq.gz > XXX/downsampled_sample_hifi.tsv
```

Use a new TSV output path for each group. The commands produce statistics from the supplied files and apply no filtering or resampling. Each Illumina pair contributes two reads. Pool all input reads before calculating N50, N90 or quality statistics.

### Result columns

| Metric | SeqKit column | Calculation |
| --- | --- | --- |
| Reads | `num_seqs` | Count of FASTQ records. |
| Total bases | `sum_len` | Sum of read lengths in bp. Divide by 1,000,000,000 for Gb. |
| N50 (bp) | `N50` | Sort reads by length from longest to shortest and report the length of the read where cumulative bases first reach 50% of the total. |
| N90 (bp) | `N90` | Apply the same calculation at 90% of total bases. The `-N 90` option adds this column. |
| Longest read (bp) | `max_len` | Maximum read length. |
| GC (%) | `GC(%)` | Percentage of all bases that are G or C. |
| Q20 (%) | `Q20(%)` | Percentage of bases with Phred score ≥20. |
| Q30 (%) | `Q30(%)` | Percentage of bases with Phred score ≥30. |
| Average quality | `AvgQual` | Convert each base's Phred score to error probability, average those probabilities over all bases, then convert the mean back to Phred using `-10 × log10(mean error probability)`. |

For a nonempty group, confirm that `min_len ≤ N90 ≤ N50 ≤ max_len`. To validate the input scope against an earlier count, compare the pooled `num_seqs` and `sum_len` with the sums of the same fields across the source files.

[SeqKit stats documentation](https://bioinf.shenwei.me/seqkit/usage/#stats) · [SeqKit quality and GC calculations](https://github.com/shenwei356/seqkit/blob/v2.11.0/seqkit/cmd/stat.go)
