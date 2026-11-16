# 06. MAG quality assessment and RNA annotation

The Workbench workflow ran these steps for each sample and strategy independently: NGS-only, ONT-only, HiFi-only, ONT-polished, HiFi-polished and NGS–ONT hybrid. Input MAGs were the raw MetaBinner bins in `.fna` format. CheckM2 and both RNA tools evaluated every raw bin. The workflow retained the raw bins and did not dereplicate across strategies.

## CheckM2 quality assessment

Workbench used CheckM2 1.1.0 with the `uniref100.KO.1.dmnd` reference database. The recorded database's base64 MD5 was `TbWEVcUf8VgEtrjhn/mfWA==`. Use the same database for all groups and place the following outputs in a new directory.

```bash
checkm2 predict --threads 8 --input XXX/raw_MAGs --extension fna \
    --database_path XXX/uniref100.KO.1.dmnd --tmpdir XXX/tmp \
    --output-directory XXX/checkm2_output
```

Check that the `Name` values in `quality_report.tsv` match the input MAG IDs exactly. Apply the following thresholds to the `Completeness` and `Contamination` percentages:

| FASTA collection | Completeness | Contamination |
| --- | --- | --- |
| Medium quality or better | ≥50% | ≤10% |
| High quality | ≥90% | <5% |

High-quality MAGs also appear in the medium-or-better set. These project thresholds use only the CheckM2 completeness and contamination values. Record both pass flags per MAG and retain the native report.

## Domain-specific rRNA and tRNA annotation

The historical binning worker wrote `MAG_domains.json` from CheckM bacterial and archaeal marker-set coverage on final bin membership. It selected the domain with stronger marker evidence and rejected missing or tied evidence. The standalone wrapper requires that validated JSON as input. Each entry contains `MAG_ID` (the FASTA basename without `.fna`) and `domain` (`bacteria` or `archaea`). The standalone binning wrapper in `05/` exports FASTAs and membership. Domain assignments must come from the historical output or an equivalent validated domain assignment step.

Use Barrnap 0.9 and tRNAscan-SE 2.0.12. The Workbench processes ran one thread per MAG and could process up to eight MAGs per group in parallel for each RNA tool.

```bash
# Bacterial MAG
barrnap -k bac --threads 1 --lencutoff 0.8 --reject 0.25 \
    XXX/bacterial_MAG.fna > XXX/bacterial_rrna.gff
tRNAscan-SE -B --thread 1 -o XXX/bacterial_trnas.tsv \
    -m XXX/bacterial_trnas_summary.txt XXX/bacterial_MAG.fna

# Archaeal MAG
barrnap -k arc --threads 1 --lencutoff 0.8 --reject 0.25 \
    XXX/archaeal_MAG.fna > XXX/archaeal_rrna.gff
tRNAscan-SE -A --thread 1 -o XXX/archaeal_trnas.tsv \
    -m XXX/archaeal_trnas_summary.txt XXX/archaeal_MAG.fna
```

Barrnap writes GFF records. The historical summary counted 5S, 16S, 23S and all rRNA records. A 16S call was counted as partial if its attributes contained `partial` or `aligned only`, and as complete otherwise. The tRNAscan-SE summary counted all calls and split them into pseudogene and non-pseudogene calls using the output annotation. The per-MAG outputs also retain the native tRNA table and run summary.

## Run the standalone wrapper

Make the three tools available on PATH. Their original Workbench tasks used separate containers. Replace every `XXX` with the relevant local directory and supply a fresh output directory.

```bash
python3 06/mag_qc.py XXX/raw_MAGs XXX/MAG_domains.json \
    XXX/uniref100.KO.1.dmnd XXX/mag_qc
```

The wrapper writes `quality_report.tsv`, `MAG_quality.csv`, `MQ_or_better_MAGs.tar.gz`, `HQ_MAGs.tar.gz`, `rrna_counts.csv`, and `trna_counts.csv`. It also keeps the native GFF, tRNA reports and logs under `annotations/` for each MAG. Workbench joined the quality and RNA counts with assembly metrics by MAG ID to produce `MAG_metrics.csv` for each group.
