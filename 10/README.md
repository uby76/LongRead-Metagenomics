# 10. Contig binning with COMEBin and SemiBin2

Each of 7 samples has 6 independent assemblies. Both binners use the same contigs strictly longer than 1,000 bp and the original sorted BAM of same-sample Illumina reads from module 05. The manifest points to the completed profile attempts. Preparation extracts only the BAM, verifies the FASTA checksum, checks reference names, lengths, order and count, and indexes the BAM once per assembly.

COMEBin 1.1.0 uses CPU training, seed 42 and two Leiden workers. The training batch size is capped at the largest power of two that fits the input, up to 1,024. SemiBin2 2.6.1 uses its bundled global model, seed 42, a minimum contig length of 1,001 bp, and the short-read algorithm for NGS-only assemblies. It uses the long-read algorithm for ONT, HiFi, polished and hybrid assemblies. Each tool retains its native bin reporting criteria.

The Workbench Nextflow entry point is `scripts/main.nf`. Supply a JSON params file with `manifest` (the URI of data/groups.tsv), `exporter` (the URI of scripts/export_bins.py), and `result` (a fresh output prefix). A recovery run may use a subset manifest and a fresh result prefix. The workflow checks the total number of successful binner outputs against the input manifest.

Results are published under `result/<sample>/<strategy>/<method>/`. The archive raw_MAGs.tar.gz contains one FASTA per raw bin. The table membership.tsv.gz includes every input contig, with an empty bin field for unassigned contigs. The table summary.tsv gives input contigs, assigned contigs and raw bin counts. Existing MetaBinner quality estimates do not apply to these new bins.

[COMEBin](https://github.com/ziyewang/COMEBin/tree/1.1.0) · [SemiBin2](https://semibin.readthedocs.io/en/latest/usage/)
