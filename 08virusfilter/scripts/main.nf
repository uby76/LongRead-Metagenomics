nextflow.enable.dsl=2
import groovy.json.JsonSlurperClassic

process FILTER_QC {
    tag { "${sample}/${strategy}" }
    cpus 8
    memory '52 GB'
    machineType 'n2-highmem-8'
    disk '140 GB'
    time '48h'
    maxForks 8
    container params.image
    publishDir "${params.result}/${sample}/${strategy}", mode: 'copy', overwrite: false
    input:
    tuple val(sample), val(strategy), val(expected), path(summary), path(fasta)
    path helper
    path database
    output:
    path 'selected.fna'
    path 'quality'
    script:
    """
    set -eu
    python '${helper}' '${summary}' '${fasta}' '${expected}'
    checkv end_to_end selected.fna quality -t ${task.cpus} -d '${database}'
    python -c "import csv; ids=[r['contig_id'] for r in csv.DictReader(open('quality/quality_summary.tsv'), delimiter=chr(9))]; src={l[1:].split()[0] for l in open('selected.fna') if l.startswith('>')}; assert len(ids)==len(src) and set(ids)==src"
    """
}

workflow {
    def rows = new JsonSlurperClassic().parseText(file(params.manifest).text)
    def jobs = rows.collect { r ->
        def base = "${params.groups}/${r.sample}/${r.strategy}/genomad_output/assembly_summary/assembly_virus"
        tuple(r.sample, r.strategy, r.expected, file(base+'_summary.tsv',checkIfExists:true), file(base+'.fna',checkIfExists:true))
    }
    FILTER_QC(Channel.fromList(jobs),Channel.value(file(params.helper,checkIfExists:true)),Channel.value(file(params.database,checkIfExists:true)))
}
