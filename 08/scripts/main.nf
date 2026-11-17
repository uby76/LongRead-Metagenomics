nextflow.enable.dsl = 2
import groovy.json.JsonSlurperClassic

process GENOMAD {
    tag { "${sample}/${strategy}" }
    cpus { size > 300000000 ? 16 : 8 }
    memory { size > 300000000 ? '112 GB' : '52 GB' }
    machineType { size > 300000000 ? 'n2-highmem-16' : 'n2-highmem-8' }
    disk '200 GB'
    time '48h'
    maxForks 8
    container params.image
    publishDir "${params.result}/groups/${sample}/${strategy}", mode: 'copy', overwrite: false, failOnError: true

    input:
    tuple val(sample), val(strategy), val(size), val(md5), path(assembly, stageAs: 'assembly.fasta.gz')
    path db, stageAs: 'db/*'
    path helper

    output:
    path 'genomad_output'
    path 'genomad.log'
    path "${sample}__${strategy}.counts.tsv", emit: counts

    script:
    """
    python '${helper}' verify-input assembly.fasta.gz '${md5}'
    python '${helper}' prepare db genomad_db
    genomad --version > genomad.log
    genomad end-to-end --threads ${task.cpus} --splits 8 assembly.fasta.gz genomad_output genomad_db >> genomad.log 2>&1
    python '${helper}' '${sample}' '${strategy}' genomad_output/assembly_summary > '${sample}__${strategy}.counts.tsv'
    """
}

process COLLECT {
    cpus 1
    memory '3 GB'
    disk '10 GB'
    time '1h'
    container params.image
    publishDir params.result, mode: 'copy', overwrite: false, failOnError: true

    input:
    path counts
    path manifest
    path helper

    output:
    path 'genomad_group_summary.tsv'

    script:
    """
    python '${helper}' merge '${manifest}' ./*.counts.tsv > genomad_group_summary.tsv
    """
}

workflow {
    def manifest = file(params.manifest, checkIfExists: true)
    def rows = new JsonSlurperClassic().parseText(manifest.text)
    assert rows.size() == 42
    assert rows.collect { [it.sample, it.strategy] }.toSet().size() == 42
    def jobs = rows.collect { row ->
        assert row.sample ==~ /[A-Za-z0-9_.-]+/
        assert row.strategy ==~ /[A-Za-z0-9_.-]+/
        tuple(row.sample, row.strategy, row.bytes, row.md5Hash, file(row.uri, checkIfExists: true))
    }
    def helper = Channel.value(file(params.helper, checkIfExists: true))
    GENOMAD(Channel.fromList(jobs), Channel.value(params.database.collect { file(it, checkIfExists: true) }), helper)
    COLLECT(GENOMAD.out.counts.collect(), Channel.value(manifest), helper)
}
