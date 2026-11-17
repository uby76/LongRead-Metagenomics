nextflow.enable.dsl = 2
import groovy.json.JsonSlurperClassic

process RGI_CONTIGS {
    tag { "${sample}/${strategy}" }
    cpus 8
    memory { bytes > 500000000 ? '112 GB' : '52 GB' }
    machineType { bytes > 500000000 ? 'n2-highmem-16' : 'n2-highmem-8' }
    disk '160 GB'
    time '48h'
    maxForks 16
    container params.image
    publishDir "${params.result}/groups/${sample}/${strategy}", mode: 'copy', overwrite: false, failOnError: true

    input:
    tuple val(sample), val(strategy), val(bytes), path(contigs)
    path card
    path runner

    output:
    path 'CARD_RGI.txt'
    path 'CARD_RGI.json'
    path 'rgi.log'
    path "${sample}__${strategy}.summary.tsv", emit: summary

    script:
    """
    echo 'dee4dcdb0d9c7f79107452d64211d816d1eab55289ddf8dc5f1e99ddfdc5e111  ${card}' | sha256sum -c -
    mkdir -p rgi_work
    cd rgi_work
    gzip -cd '../${contigs}' > input.fasta
    {
        rgi main --version
        rgi load --card_json '../${card}' --local
        rgi main -i input.fasta -o CARD_RGI -t contig -a DIAMOND -n ${task.cpus} --data wgs --local
    } > ../rgi.log 2>&1
    python '../${runner}' '${sample}' '${strategy}' CARD_RGI > '../${sample}__${strategy}.summary.tsv'
    cp CARD_RGI.txt CARD_RGI.json ..
    """
}

process SUMMARIZE {
    cpus 1
    memory '3 GB'
    disk '10 GB'
    time '1h'
    container params.image
    publishDir params.result, mode: 'copy', overwrite: false, failOnError: true

    input:
    path tables
    path manifest
    path runner

    output:
    path 'CARD_RGI_group_summary.tsv'

    script:
    """
    python '${runner}' merge '${manifest}' ./*.summary.tsv > CARD_RGI_group_summary.tsv
    """
}

workflow {
    def manifest = file(params.manifest, checkIfExists: true)
    def rows = new JsonSlurperClassic().parseText(manifest.text)
    assert rows.size() == 42
    assert rows.collect { [it.sample, it.strategy] }.toSet().size() == 42
    def inputs = rows.collect { row -> tuple(row.sample, row.strategy, row.bytes, file(row.uri, checkIfExists: true)) }
    def runner = Channel.value(file(params.runner, checkIfExists: true))
    RGI_CONTIGS(Channel.fromList(inputs), Channel.value(file(params.card, checkIfExists: true)), runner)
    SUMMARIZE(RGI_CONTIGS.out.summary.collect(), Channel.value(manifest), runner)
}
