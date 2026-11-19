# 09. Cross-strategy viral diversity and assembly quality

Run from the repository root with Python 3, requests, gcloud and Workbench storage access. Choose a new snapshot directory for each refresh.

    python3 09/scripts/preflight.py 09/data/snapshot_NAME

The snapshot records publication status for every expected group, preserves the sample, assembly strategy, original contig, geNomad virus ID and coordinates, and joins existing geNomad taxonomy and CheckV quality by exact sequence ID. Sequence lengths and expected counts must agree. Counts use the number of selected records in each group as denominator. A pending group has missing data. Comparisons across strategies should use the same samples in each strategy and report missing completeness separately.

Once all 42 groups are published, run on a VM with BLAST+, sort, gcloud, scratch space and access to the original objects. The command copies each input to temporary storage and preserves only one representative FASTA and the membership table.

    python3 09/scripts/cluster.py 09/data/snapshot_NAME 09/data/votu_NAME --threads 8

The implementation uses all-versus-all megablast. For each shorter query/longer subject pair it greedily retains HSPs by bit score, dropping any HSP that overlaps retained coordinates on either sequence. ANI is summed identical positions divided by summed alignment lengths. Coverage uses covered query bases divided by shorter sequence length. The inclusive cutoffs are 95% ANI and 85% coverage. Longest-first greedy representative clustering requires a direct qualifying alignment to the representative. This conservative overlap handling can exclude borderline pairs. No full production clustering run has yet validated this implementation.

The first snapshot is stored in the longshort project under jianwei_result/09/data/snapshot_20261118/. CheckV and geNomad results are reused. The sequence table contains the original geNomad lineage for each member. Monodnaviria alone does not establish a bacterial host. Host linkage, lifestyle, functional annotation and vConTACT3 outputs are not present in this release.

## Published-group reports

Run scripts/group_report.py with the preflight snapshot path, a new report directory, and the original --rgi-manifest and --genomad-manifest JSON files. It checks sample, strategy, assembly URI and size before joining. The original RGI manifest lacks generation and checksum fields, so this comparison does not establish historical byte identity.

quality/by_group.tsv reports all selected records, missing completeness values and median CheckV completeness among available values. quality/six_strategy_matched.tsv retains the same samples for each strategy. taxonomy/genomad.tsv counts native geNomad lineages. arg/candidate_counts.tsv counts Perfect/Strict RGI ORFs fully contained within a selected viral region in the same sample and strategy. Provirus boundaries use the original host-contig coordinates, with both endpoints included. Boundary-crossing ORFs are excluded. Counts are per assembly and repeated recovery across strategies contributes multiple records. These are sequence-based candidates; gene function and viral carriage have not been experimentally established.

The second interim report is under jianwei_result/09/group_report_20261118_02/. Its input_status.tsv freezes the group denominator. Each input-status row points to the upstream FASTA. Completeness and taxonomy results can be updated by running the same scripts in a fresh directory after publication completes.

The managed Workbench entry point is scripts/votu.nf. It stages the 42 selected FASTAs in input_status.tsv order. Outputs in 09/vOTU include membership.tsv (source sample, strategy, contig, geNomad coordinates and member ID), representatives.fna, genomad_taxonomy.tsv (representative lineage and number of distinct member lineages), and versions.txt. Source lineage calls are reused.
