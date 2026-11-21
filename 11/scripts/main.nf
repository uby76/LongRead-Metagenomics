nextflow.enable.dsl=2

process DEREPLICATE {
    cpus 32
    memory '190 GB'
    machineType 'n2-highmem-32'
    disk '300 GB'
    time '72h'
    container params.image
    publishDir "${params.result}/dereplication", mode:'copy', overwrite:false, failOnError:true
    input:
    path bundles, stageAs:'inputs/group??/MQ.tar.gz'
    path metrics
    path worker
    output:
    path 'membership.tsv'
    path 'representatives.tar.gz', emit: representatives
    path 'skani_version.txt'
    script:
    """
    set -euo pipefail
    python '$worker' unpack '$metrics'
    skani --version > skani_version.txt
    skani triangle --sparse --min-af 49 -t ${task.cpus} genomes/*.fna -o pairs.tsv
    python '$worker' cluster '$metrics'
    """
}

process REFERENCE {
    cpus 4
    memory '12 GB'
    disk '140 GB'
    time '48h'
    errorStrategy 'retry'
    maxRetries 1
    container params.image
    output:
    path 'gtdbtk_r232_data.tar.gz'
    script:
    """
    set -eu
    python -c 'import urllib.request; urllib.request.urlretrieve("https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz","gtdbtk_r232_data.tar.gz")'
    echo '25a59e0352b1fd150c589f56559767d4  gtdbtk_r232_data.tar.gz' | md5sum -c -
    """
}

process TAXONOMY_AND_TREES {
    cpus 32
    memory '230 GB'
    machineType 'n2-highmem-32'
    disk '1000 GB'
    time '168h'
    container params.image
    publishDir "${params.result}/taxonomy_phylogeny", mode:'copy', overwrite:false, failOnError:true, saveAs:{ it.tokenize('/')[-1] }
    input:
    path representatives
    path db
    path worker
    output:
    path 'results/*'
    script:
    """
    set -euo pipefail
    mkdir reps database results
    tar -xzf '$representatives' -C reps
    tar -xzf '$db' -C database
    metadata=\$(find database -type f -path '*/metadata/metadata.txt')
    test -f "\$metadata"
    export GTDBTK_DATA_PATH=\$(realpath "\$(dirname "\$(dirname "\$metadata")")")
    gtdbtk --version > results/gtdbtk_version.txt
    gtdbtk classify_wf --genome_dir reps -x fna --out_dir classify --cpus ${task.cpus} --pplacer_cpus 4 --place_species --keep_intermediates
    python '$worker' finish
    for f in classify/classify/*.summary.tsv classify/align/*.user_msa.fasta.gz classify/align/*.filtered.tsv classify/identify/*.failed_genomes.tsv; do
        if test -f "\$f"; then cp -L "\$f" results/; fi
    done
    """
}

workflow {
    def archives = file(params.archives, checkIfExists:true).readLines().collect { file(it,checkIfExists:true) }
    assert archives.size()==42 && archives.toSet().size()==42
    def worker = file(params.worker,checkIfExists:true)
    DEREPLICATE(Channel.value(archives),file(params.metrics,checkIfExists:true),worker)
    REFERENCE()
    TAXONOMY_AND_TREES(DEREPLICATE.out.representatives,REFERENCE.out.first(),worker)
}
