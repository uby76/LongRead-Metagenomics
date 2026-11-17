# 07. CARD/RGI annotation of assembled contigs

This directory packages the CARD/RGI Nextflow workflow used by the longshort project. It annotates assembled contigs for seven samples (CNTL, FF, FFPS, KB, KBPS, YGC and YML) and six strategies (NGS-only, ONT-only, HiFi-only, ONT-polished, HiFi-polished and NGS–ONT hybrid): 42 independent groups. Inputs are gzip-compressed contig FASTAs, not the individual MAG bins or raw reads.

## Files

- `scripts/main.nf`: original two-process Nextflow DSL2 workflow, unchanged.
- `scripts/validate.py`: original output validation and group-summary script, unchanged.
- `scripts/nextflow.config`: original execution settings with an optional local Docker profile.
- `data/manifest.example.json`: all 42 sample/strategy entries, with portable input-path placeholders. The `bytes` values are retained from the historical manifest.
- `data/params.example.json`: input, output and pinned-container settings.
- `plots/`: reserved for final PDF figures; this workflow does not generate plots.

Source: the longshort `CARD_RGI_contigs/code` workflow registered on 16 November 2026, locally retained under `MAGs/MAGs6/03CARD_RGI/code`. Cloud account details, submission records, input sequences, database files and analysis outputs are not bundled.

## Reference database and method

Use CARD **4.0.1** `card.json`, with SHA-256:

```text
dee4dcdb0d9c7f79107452d64211d816d1eab55289ddf8dc5f1e99ddfdc5e111
```

This checksum was verified against the locally retained database. The workflow checks it in every annotation task. Obtain that release from CARD and supply the extracted `card.json`; a newer database will intentionally fail the check. The original RGI container is pinned in the parameter template:

```text
quay.io/biocontainers/rgi@sha256:164f4a58042146d144af35bc92589b47e095e40b4588836363e4b6025ac014b4
```

The RGI version is recorded at runtime in each `rgi.log`. The original commands are:

```bash
rgi main --version
rgi load --card_json ../card.json --local
rgi main -i input.fasta -o CARD_RGI -t contig -a DIAMOND -n 8 --data wgs --local
```

Each task uses its own local database/work directory. No `--include_loose` or `--low_quality` option is added. Validation accepts only `Perfect` and `Strict` calls and reports their counts; it does not apply an extra identity or coverage threshold.

## Run

Install Nextflow and a compatible container runtime. From the repository root, edit the two JSON templates, replacing every `XXX`. Keep the complete 42-entry manifest for the historical design. Each entry supplies `sample`, `strategy`, `uri` (a local path or accessible cloud URI), and `bytes` (compressed input file size). Update `bytes` if the input file changes. Use unique sample/strategy pairs and simple filesystem-safe identifiers.

For a local Docker run:

```bash
nextflow run 07/scripts/main.nf \
    -c 07/scripts/nextflow.config \
    -params-file 07/data/params.example.json \
    -profile docker
```

Use a fresh result directory. The workflow publishes with `overwrite: false`, retains its work directory and allows one retry. To resume the same run with the same inputs and work directory, add `-resume`.

For Google Batch, supply your own executor/project/location configuration with an additional `-c` file, use accessible cloud paths and omit `-profile docker`. The original managed Workbench executor configuration and credentials are not part of this package.

The historical workflow requests eight CPUs per annotation group, 52 GB RAM for inputs up to 500,000,000 compressed bytes or 112 GB above that size, 160 GB disk and 48 hours. Google Batch machine types are `n2-highmem-8` and `n2-highmem-16`, respectively. Up to 16 annotation tasks may run concurrently. The summary process requests one CPU, 3 GB RAM, 10 GB disk and one hour. Adjust concurrency for available infrastructure without changing annotation arguments.

## Outputs and validation

Each `groups/<sample>/<strategy>/` directory contains `CARD_RGI.txt`, `CARD_RGI.json`, `rgi.log` and a group summary TSV. The result root contains `CARD_RGI_group_summary.tsv`, with `sample_ID`, `sequencing_strategy`, `perfect_hits`, `strict_hits` and `total_hits`.

The original validator checks JSON readability, the `Cut_Off` column, allowed hit categories and complete coverage of the 42 expected groups. Counts are annotation rows, not deduplicated ARG families, abundance estimates or copies per cell. The validator uses Python assertions: do not run it with `python -O` or `PYTHONOPTIMIZE` enabled.

Packaging checks cover Python syntax, configuration/manifest consistency and synthetic validation cases. The full Nextflow/RGI analysis was not rerun during packaging; raw contigs and the reference database are required to execute it.
