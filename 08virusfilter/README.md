# geNomad virus selection and CheckV assessment

Across 42 sample/strategy groups the geNomad lineage contains the exact token Monodnaviria in 1,489 predictions. The spelling monodnaviridae has zero matches (case-insensitive). The audit counts 35,436 predictions longer than 10,000 bp and another 904 Monodnaviria predictions from 1,001 through 10,000 bp. Counts represent records across assemblies.

The workflow reads assembly_virus_summary.tsv and assembly_virus.fna from each geNomad group. It selects from the native viral FASTA, preserving geNomad provirus boundaries, and verifies sequence IDs, lengths and expected counts. CheckV 1.0.3 uses database v1.5. Each group retains selected.fna and the CheckV quality directory, including all quality grades in quality_summary.tsv.

Set the paths in data/params.example.json. From the repository root run Nextflow with a configured container executor, for example:

    nextflow run 08virusfilter/scripts/main.nf -params-file 08virusfilter/data/params.example.json -with-docker

scripts/audit.py accepts native virus summary TSVs and reports length and lineage counts. The archived audit uses sample__strategy.tsv filenames.
