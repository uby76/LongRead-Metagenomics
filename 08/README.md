# 08. geNomad annotation of assembled contigs

This directory packages the geNomad workflow submitted for the longshort project on 17 November 2026: seven samples (CNTL, FF, FFPS, KB, KBPS, YGC and YML) × six assembly strategies (NGS-only, ONT-only, HiFi-only, ONT-polished, HiFi-polished and NGS–ONT hybrid), totaling 42 independent groups.

## Scope and method

Use the complete gzip-compressed contig FASTAs from each assembly's `build/contigs.fasta.gz`, before the length filtering used for binning. No additional input length threshold is applied. Existing assembler output limits and polishing still affect the available sequences.

```bash
genomad end-to-end --threads 8 --splits 8 assembly.fasta.gz genomad_output genomad_db
```

The workflow uses 16 threads for compressed inputs larger than 300,000,000 bytes and eight otherwise. It retains geNomad's default classification filters, including its native length-dependent criteria; neither `--relaxed` nor `--conservative` is added. Provirus identification remains enabled. All geNomad output directories, including intermediate files, predicted virus/plasmid sequences, taxonomy, provirus coordinates and logs, are retained.

**This workflow runs geNomad only.** The proposed downstream >10 kbp rule, >1 kbp taxonomy exception and CheckV assessment are deferred and are not implemented here. No plots are produced.

## Files and dependencies

- `scripts/main.nf`: submitted Nextflow DSL2 workflow, unchanged.
- `scripts/results.py`: submitted database-link preparation, input-checksum validation, output validation and group-count aggregation, unchanged.
- `scripts/nextflow.config`: submitted execution settings plus an optional local Docker profile.
- `data/manifest.example.json`: 42 input entries with portable paths and historical compressed sizes/base64 MD5 values.
- `data/params.example.json`: input, database, output and container settings.
- `data/reference_manifest.json`: file sizes and checksums of the database used for the submitted analysis.
- `plots/`: reserved for final PDF figures.

Use Nextflow and a compatible container runtime. The original container is pinned by digest:

```text
quay.io/biocontainers/genomad@sha256:ef01654b515ca5d9489eee57d36b1571e8b73396d267198c014ba7744a0a87b5
```

The executable version is logged in each group's `genomad.log`. Use geNomad database **1.7**. The helper checks `version.txt` and reconstructs mini-database symbolic links when cloud storage has represented those links as short text files. Keep all database files listed in the parameter template. The reference manifest describes the archived cloud representation: link-placeholder checksums can differ from those of extracted native symbolic links. Database checksums are recorded for provenance; they are not automatically enforced by the workflow. Do not enable Python optimization (`-O` or `PYTHONOPTIMIZE`), which disables the database-version assertion.

## Run

From the repository root, replace every `XXX` in both example JSON files with the appropriate path. Input paths may be local files or accessible cloud URIs. Use a fresh result directory. The workflow intentionally requires the complete 42-group design and unique sample/strategy pairs.

If supplying different input files, update both `bytes` and `md5Hash`; the latter must be the **base64-encoded binary MD5**, not hexadecimal MD5. Obtain it for each compressed FASTA with:

```bash
openssl dgst -md5 -binary XXX/contigs.fasta.gz | openssl base64 -A
```

For a local Docker run:

```bash
nextflow run 08/scripts/main.nf \
    -c 08/scripts/nextflow.config \
    -params-file 08/data/params.example.json \
    -profile docker
```

For Google Batch, omit `-profile docker` and supply your own executor/project/location configuration with an additional `-c` file. The managed Workbench configuration, account details, raw sequences, reference database binaries and submission records are not bundled. This GitHub upload does not submit another analysis.

The workflow allows up to eight concurrent annotation tasks. Tasks request 52 GB RAM/eight CPUs or 112 GB RAM/16 CPUs for compressed inputs >300 MB, 200 GB disk and 48 h. The corresponding Google Batch machine types are `n2-highmem-8` and `n2-highmem-16`. The aggregation task requests one CPU, 3 GB RAM, 10 GB disk and one hour. Tasks allow one retry, published files are not overwritten, and the work directory is retained. Add `-resume` to resume the same run with unchanged inputs and work directory. Set concurrency to suit available infrastructure.

## Outputs and checks

Each `groups/<sample>/<strategy>/` directory receives `genomad_output/`, `genomad.log` and a group counts TSV. The result root receives `genomad_group_summary.tsv` after every group completes. Columns are `sample_ID`, `sequencing_strategy`, `viral_predictions`, `provirus_predictions` and `plasmid_predictions`. Provirus predictions are a subset of viral predictions, not an additional category to sum.

Before annotation, the helper checks the compressed input's MD5 against the manifest. After annotation it checks unique sequence IDs and exact FASTA/summary length agreement for both virus and plasmid predictions. Aggregation requires all expected groups with no duplicates.

Packaging checks cover Python syntax, script identity with the submitted workflow and JSON consistency. Synthetic helper checks performed before submission covered sequence-length mismatches, empty plasmid output, provirus counts, duplicate-group rejection and input MD5 validation. These checks do not establish completion of the full cloud analysis.

See the [geNomad quickstart](https://portal.nersc.gov/genomad/quickstart.html) and [provirus identification documentation](https://portal.nersc.gov/genomad/provirus_identification.html).
