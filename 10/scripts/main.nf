nextflow.enable.dsl=2
params.prepared = null

process PREPARE {
    tag "${sample}/${strategy}"
    cpus 4
    machineType 'n2-standard-4'
    disk '500 GB'
    errorStrategy 'ignore'
    maxForks 4
    container 'quay.io/biocontainers/metabinner@sha256:e700129938b2f12473e00563474d56fdb6173bf8e93f401a8e75ab34a91616a2'
    input:
    tuple val(sample), val(strategy), val(count), val(digest), path(fasta), path(archive)
    output:
    tuple val(sample), val(strategy), val(count), path('contigs.fa'), path('reads.bam'), path('reads.bam.bai')
    script:
    """
    set -euo pipefail
    echo '$digest  $fasta' | sha256sum -c -
    gzip -dc '$fasta' > contigs.fa
    python3 - '$archive' <<'PYTHON'
import shutil
import sys
import tarfile
from pathlib import PurePosixPath

found = False
with tarfile.open(sys.argv[1], 'r|gz', bufsize=1024 * 1024) as archive:
    for entry in archive:
        if entry.isfile() and PurePosixPath(entry.name).match('coverage_output/work_files/*.bam'):
            if found:
                raise ValueError('Multiple archived BAMs')
            with open('reads.bam', 'wb') as output:
                shutil.copyfileobj(archive.extractfile(entry), output, 1024 * 1024)
            found = True
if not found:
    raise ValueError('Missing archived BAM')
PYTHON
    samtools quickcheck -v reads.bam
    samtools faidx contigs.fa
    samtools view -H reads.bam | awk -F '\\t' '\$1=="@SQ"{for(i=2;i<=NF;i++){if(\$i~/^SN:/)n=substr(\$i,4);if(\$i~/^LN:/)l=substr(\$i,4)}print n "\\t" l}' > bam.refs
    cut -f1,2 contigs.fa.fai > fasta.refs
    cmp fasta.refs bam.refs
    test "\$(wc -l < fasta.refs)" -eq '$count'
    samtools index -@ 3 reads.bam
    """
}

process BIN {
    tag "${method}/${sample}/${strategy}"
    cpus 16
    machineType { method == 'SemiBin2' && strategy == 'NGS-only' && count.toInteger() > 250000 ? 'n2-highmem-96' : count.toInteger() > 250000 ? 'n2-highmem-32' : 'n2-highmem-16' }
    disk '400 GB'
    maxForks 6
    time '168h'
    errorStrategy 'ignore'
    container { method == 'COMEBin' ? 'quay.io/biocontainers/comebin@sha256:32ce04758ca5d19450809d2423f2b4a28e50687c7a9f565affa759ba247f9c6f' : 'quay.io/biocontainers/semibin@sha256:9f4c3b0c58483053b93ded62bb73dcf47e47c25a7c9e557358103e3e69a74142' }
    publishDir { "${params.result}/${sample}/${strategy}/${method}" }, mode: 'copy', overwrite: false, failOnError: true
    input:
    tuple val(sample), val(strategy), val(count), path(fasta), path(bam), path(bai), val(method)
    path exporter
    output:
    path 'raw_MAGs.tar.gz'
    path 'membership.tsv.gz'
    path 'summary.tsv', emit: summaries
    script:
    def mode = strategy == 'NGS-only' ? 'short_read' : 'long_read'
    def batch = Math.min(1024, Integer.highestOneBit(count.toInteger()))
    def command = method == 'COMEBin' ?
        "mkdir align && ln -s ../reads.bam ../reads.bam.bai align/ && run_comebin.sh -a contigs.fa -p align -o run -t ${task.cpus} -w 2 -d cpu -s 42 -b ${batch}" :
        "SemiBin2 single_easy_bin -i contigs.fa -b reads.bam -o run --environment global --sequencing-type ${mode} --min-len 1001 --engine cpu --random-seed 42 -t ${task.cpus} --compression none"
    def bins = method == 'COMEBin' ? 'run/comebin_res/comebin_res_bins' : 'run/output_bins'
    """
    set -euo pipefail
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=1
    if [ '${method}' = 'COMEBin' ]; then
        mkdir compat
        cat > compat/realpath <<'PYTHON'
#!/usr/bin/env python3
import argparse
import os
parser = argparse.ArgumentParser()
parser.add_argument('-e', action='store_true')
parser.add_argument('-m', action='store_true')
parser.add_argument('paths', nargs='+')
args = parser.parse_args()
for path in args.paths:
    if args.e and not os.path.exists(path):
        raise SystemExit(1)
    print(os.path.realpath(path))
PYTHON
        chmod +x compat/realpath
        export PATH="\$PWD/compat:\$PATH"
    fi
    ${command}
    python "${exporter}" "${sample}" "${strategy}" "${method}" contigs.fa "${bins}"
    """
}
workflow {
    def rows = Channel.fromPath(params.manifest, checkIfExists:true).splitCsv(header:true, sep:'\t')
    def inputs = rows.map { r -> tuple(r.sample, r.strategy, r.contigs, r.sha256, file(r.fasta), file(r.archive)) }
    def expected = rows.count().map { it * 2 }
    PREPARE(inputs)
    def fresh = PREPARE.out.combine(Channel.of('COMEBin','SemiBin2'))
    def reused = params.prepared ? Channel.fromPath(params.prepared, checkIfExists:true).splitCsv(header:true, sep:'\t').map { r ->
        tuple(r.sample, r.strategy, r.contigs, file(r.prefix+'/contigs.fa'), file(r.prefix+'/reads.bam'), file(r.prefix+'/reads.bam.bai'), r.method)
    } : Channel.empty()
    expected = expected.combine(reused.count()).map { a,b -> a+b }
    BIN(fresh.mix(reused), file(params.exporter, checkIfExists:true))
    BIN.out.summaries.count().combine(expected).subscribe { done, total ->
        if (done != total || !total) error("Completed ${done} of ${total} binning analyses")
    }
}
