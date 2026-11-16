# 02. Contig assembly

These commands reproduce the assembler settings recorded for the primary single-platform assemblies of the seven matched samples. Each sample and platform was assembled separately using its downsampled input. Replace each `XXX` with the appropriate directory. Run the wrapper commands from the repository root and use a new assembler output directory for each sample.

## Illumina paired-end assembly

MEGAHIT v1.2.9 used the downsampled, deduplicated R1 and R2 FASTQ files for the same sample. The recorded command used the `meta-large` preset, 16 threads, a memory limit of 412,316,860,416 bytes (384 GiB), and a minimum contig length of 300 bp.

```bash
megahit -1 XXX/sample_i_ds_1.dedup.fastq.gz -2 XXX/sample_i_ds_2.dedup.fastq.gz \
    --presets meta-large --min-contig-len 300 -t 16 -m 412316860416 \
    -o XXX/megahit_out
```

MEGAHIT reported the following effective k-mer lengths in all seven assembly logs:

```text
27,37,47,57,67,77,87,97,107,117,127
```

The log recorded that the maximum k was reset to 127 for these inputs. The original command selected the preset without an explicit `--k-list`.

MEGAHIT writes the assembly to `XXX/megahit_out/final.contigs.fa` and the log to `XXX/megahit_out/log`. The historical pipeline published the gzip-compressed FASTA as `sample_i_ds.contigs.fasta.gz` and retained a sample-named MEGAHIT log. The wrapper creates `contigs.fasta.gz` alongside the native files.

```bash
bash 02/assemble_illumina.sh XXX/sample_i_ds_1.dedup.fastq.gz \
    XXX/sample_i_ds_2.dedup.fastq.gz XXX/megahit_out
```

## Oxford Nanopore and PacBio HiFi assembly

MetaMDBG v1.4 took one cleaned, downsampled FASTQ per sample. ONT used `--in-ont` and PacBio HiFi used `--in-hifi`. Both used eight threads and a minimum contig length of 300 bp.

```bash
# Oxford Nanopore
metaMDBG --out-dir XXX/sample_n_ds_mdbg_out \
    --in-ont XXX/sample_n_ds.cleaned.fastq.gz --threads 8 --min-contig-length 300

# PacBio HiFi
metaMDBG --out-dir XXX/sample_p_ds_mdbg_out \
    --in-hifi XXX/sample_p_ds.cleaned.fastq.gz --threads 8 --min-contig-length 300
```

| Setting reported in the logs | ONT | PacBio HiFi |
| --- | --- | --- |
| Minimizer length (bp) | 15 | 15 |
| Correction density | 0.025 | 0.025 |
| Assembly density | 0.005 | 0.005 |
| Read correction | On | Off |
| Homopolymer compression | Off | On |

The 15 bp value is the length of an individual minimizer. MetaMDBG builds its graph with tuples of minimizers. The original commands used the version's minimizer defaults.

MetaMDBG writes `contigs.fasta.gz` and `metaMDBG.log` inside each output directory. The historical pipeline also published the compressed FASTA under the sample name, such as `sample_n_ds.contigs.fasta.gz` or `sample_p_ds.contigs.fasta.gz`.

```bash
bash 02/assemble_long_reads.sh ont XXX/sample_n_ds.cleaned.fastq.gz \
    XXX/sample_n_ds_mdbg_out
bash 02/assemble_long_reads.sh hifi XXX/sample_p_ds.cleaned.fastq.gz \
    XXX/sample_p_ds_mdbg_out
```
