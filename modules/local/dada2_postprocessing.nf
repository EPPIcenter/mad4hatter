// Define a Nextflow module for creating a sequence table from DADA2 mergers
process DADA2_POSTPROCESSING {

    label 'process_medium'

    // Define input channels

    input:
    path seqtab_paths
    path amplicon_file
    val verbose
    val bimera_removal_method

    // Define output channels
    output:
    path("dada2.clusters.txt"), emit: dada2_clusters

    script:
    def seqtab_paths = seqtab_paths.join(' ')
    """
    Rscript ${projectDir}/bin/dada2_postprocessing.R \
        --seqtab-paths ${seqtab_paths} \
        --amplicon-file ${amplicon_file} \
        --ncores ${task.cpus} \
        --bimera-removal-method ${bimera_removal_method} \
        --log-level ${params.logLevel}
    """
}
