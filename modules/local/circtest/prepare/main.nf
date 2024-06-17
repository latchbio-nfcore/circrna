process CIRCTEST_PREPARE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biocontainers/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(circ_counts)
    tuple val(meta2), path(gene_counts)

    output:
    tuple val(meta), path('*_circs.tsv'), emit: circ_counts
    tuple val(meta), path('*_genes.tsv'), emit: gene_counts

    path "versions.yml"                 , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template 'prepare.R'
}
