# 05. MetaBinner contig binning

The Workbench workflow applied these steps independently to each sample and strategy: NGS-only, ONT-only, HiFi-only, ONT-polished, HiFi-polished, and NGS–ONT hybrid. Each group used its own final contig FASTA and the same sample's cleaned, downsampled Illumina R1/R2 files for coverage. Each group generated its own coverage and composition profiles.

## Run one group

Run from the repository root in a Linux environment with MetaBinner and its dependencies, including BWA, SAMtools and the JGI depth tool. Set the installation directory to the directory containing `run_metabinner.sh` and `scripts/`. Use a new output directory. Use paths without whitespace because the upstream coverage script expects unquoted paths.

```bash
bash 05/metabinner.sh XXX/MetaBinner XXX/contigs.fasta.gz \
    XXX/sample_i_ds_1.cleaned.fastq.gz XXX/sample_i_ds_2.cleaned.fastq.gz \
    XXX/metabinner_group sample__ONT-polished
```

Use a group ID that identifies the sample and strategy. The Workbench runs used the MetaBinner container pinned to `quay.io/biocontainers/metabinner@sha256:e700129938b2f12473e00563474d56fdb6173bf8e93f401a8e75ab34a91616a2`.

## Filter contigs and calculate profiles

The worker staged an uncompressed `source_contigs.fasta` from the selected assembly. The official filter retained contigs **strictly longer than 1,000 bp**. The wrapper verifies unique original contig IDs and checks that the filtered FASTA contains exactly the expected IDs in order.

```bash
python3 XXX/MetaBinner/scripts/Filter_tooshort.py XXX/source_contigs.fasta 1000
cp XXX/source_contigs_1000.fa XXX/contigs_filtered.fasta

python3 XXX/MetaBinner/scripts/gen_kmer.py XXX/contigs_filtered.fasta 1000 4

bash XXX/MetaBinner/scripts/gen_coverage_file.sh \
    -a XXX/contigs_filtered.fasta -o XXX/coverage_output \
    -t 16 -m 64 -l 1000 -f _1.cleaned.fastq.gz -r _2.cleaned.fastq.gz \
    XXX/sample_i_ds_1.cleaned.fastq.gz XXX/sample_i_ds_2.cleaned.fastq.gz
```

The 4-mer script generates `contigs_filtered_kmer_4_f1000.csv`. The coverage script aligns paired reads with BWA MEM, sorts the alignments with SAMtools and calculates contig depth with `jgi_summarize_bam_contig_depths`. The worker used `coverage_output/coverage_profile.tsv`. The Workbench worker preserved generated SAM files. This standalone wrapper uses the official coverage script's standard intermediate-file handling.

The wrapper checks that the filtered FASTA and both profiles contain exactly the same unique contig IDs, then writes both profiles in FASTA order. Binning uses those reordered tables.

## Run MetaBinner and export MAGs

```bash
bash XXX/MetaBinner/run_metabinner.sh -a XXX/contigs_filtered.fasta \
    -o XXX/metabinner_output -d XXX/coverage_profile.tsv \
    -k XXX/kmer_profile.csv -p XXX/MetaBinner -s large -t 8
```

The final assignment file is `metabinner_output/metabinner_res/metabinner_result.tsv`. Each nonblank line contains a contig ID and bin ID, separated by a tab. The exporter verifies that every assigned contig exists in the filtered FASTA and occurs in only one final assignment. It writes one `raw_MAGs/<group_id>__bin_<bin_id>.fna` per bin, plus `raw_MAGs.tar.gz`, `final_contig_to_MAG.csv` and `unbinned_contigs.csv`. These FASTAs contain the raw bins before CheckM2 quality selection.

[MetaBinner upstream repository](https://github.com/ziyewang/MetaBinner)
