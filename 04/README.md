# 04. Hybrid metagenome assembly with metaSPAdes

Use SPAdes 4.3.0 and replace each `XXX` with the appropriate directory. Supply synchronized Illumina mate files and a long-read FASTQ from the same sample. Use a separate output directory for each sample and strategy.

## Illumina–ONT workflow recorded in Workbench

The run used the cleaned, downsampled Illumina R1 and R2 files and the cleaned, downsampled ONT FASTQ. The worker checked that the Illumina mates were paired before assembly. The command below matches the verified YML log.

```bash
spades.py --meta --only-assembler \
    -1 XXX/sample_i_ds_1.cleaned.fastq.gz \
    -2 XXX/sample_i_ds_2.cleaned.fastq.gz \
    --nanopore XXX/sample_n_ds.cleaned.fastq.gz \
    -t 32 -m 256 --checkpoints all -o XXX/ont_hybrid
```

The command configured **32 threads** and a **256 GB** memory limit. It used SPAdes automatic k-mer selection. The YML log reported **k = 21, 33, 55** and recorded those stages. Check `spades.log` for the effective k values for each sample. The workflow used `--only-assembler` to start assembly from the previously cleaned reads.

## Illumina–PacBio HiFi example

The project records contain completed Illumina–ONT assemblies. A historical Illumina–HiFi hybrid run, its effective k values and its resource settings could not be verified. The following example adapts the ONT command to SPAdes' PacBio input option and uses the same resource settings as a starting point. Hybrid metaSPAdes is experimental, so validate the example before using it for a batch.

```bash
spades.py --meta --only-assembler \
    -1 XXX/sample_i_ds_1.cleaned.fastq.gz \
    -2 XXX/sample_i_ds_2.cleaned.fastq.gz \
    --pacbio XXX/sample_p_ds.cleaned.fastq.gz \
    -t 32 -m 256 --checkpoints all -o XXX/hifi_hybrid
```

The HiFi example also uses automatic k-mer selection. Record the effective k list from the new run's log.

## Checkpoints and thread verification

Start the run with `--checkpoints all`. Preserve the full output directory, including `params.txt`, `pipeline_state/` and the k-stage files. When restoring a snapshot on another machine, make the same inputs available at the paths stored by SPAdes. Use the same SPAdes version and run the continuation from the existing output directory by passing its path to `-o`.

```bash
OMP_NUM_THREADS=32 OMP_THREAD_LIMIT=32 OMP_DYNAMIC=FALSE \
    spades.py --continue -o XXX/ont_hybrid
```

The Workbench recovery restored a native checkpoint and applied these OpenMP settings to the saved 32-thread run. Confirm the effective thread count in the SPAdes core log line containing `Maximum # of threads to use (adjusted due to OMP capabilities)`. The same resume procedure applies to the HiFi example if a valid checkpoint exists.

## Run the wrapper

```bash
bash 04/hybrid_assembly.sh ont XXX/sample_i_ds_1.cleaned.fastq.gz \
    XXX/sample_i_ds_2.cleaned.fastq.gz XXX/sample_n_ds.cleaned.fastq.gz XXX/ont_hybrid
bash 04/hybrid_assembly.sh hifi XXX/sample_i_ds_1.cleaned.fastq.gz \
    XXX/sample_i_ds_2.cleaned.fastq.gz XXX/sample_p_ds.cleaned.fastq.gz XXX/hifi_hybrid
bash 04/hybrid_assembly.sh resume XXX/ont_hybrid
```

Use `contigs.fasta` as the assembly output and keep `spades.log` for parameters and stage history. Workbench verified a nonempty contig FASTA, published a compressed copy as `contigs.fasta.gz`, and saved assembly statistics. The wrapper verifies the native FASTA and leaves it in the SPAdes output directory.

[SPAdes command-line options](https://ablab.github.io/spades/running.html)
