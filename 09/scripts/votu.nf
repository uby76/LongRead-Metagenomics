nextflow.enable.dsl=2

process CLUSTER {
    tag '42 assemblies'
    cpus 32
    memory '115 GB'
    machineType 'n2-standard-32'
    disk '400 GB'
    time '96h'
    container params.image
    publishDir params.result, mode: 'copy', overwrite: false, failOnError: true, saveAs: { it.tokenize('/')[-1] }
    input:
    path fasta, stageAs: 'assemblies/group????/selected.fna'
    path status, stageAs: 'snapshot/input_status.tsv'
    path metadata, stageAs: 'snapshot/sequence_metadata.tsv'
    path script_file
    output:
    path 'output/membership.tsv'
    path 'output/representatives.fna'
    path 'output/genomad_taxonomy.tsv'
    path 'versions.txt'
    script:
    """
    set -eu
    export LC_ALL=C
    export TMPDIR=\$PWD
    { blastn -version; makeblastdb -version; python --version; } > versions.txt
    python '${script_file}' snapshot output --staged assemblies --threads ${task.cpus}
    """
}
workflow {
    def status = file(params.status,checkIfExists:true)
    def rows = status.readLines().drop(1).collect { it.split('\t') }
    assert rows.size()==42 && rows.every { it[3]=='published' }
    CLUSTER(Channel.value(rows.collect { file(it[4],checkIfExists:true) }),Channel.value(status),Channel.value(file(params.metadata,checkIfExists:true)),Channel.value(file(params.script_file,checkIfExists:true)))
}
