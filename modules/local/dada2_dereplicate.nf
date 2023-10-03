// Define a Nextflow module for Dereplication using DADA2
process DADA2_DEREPLICATE {

    label 'process_low'

    // Define input channels
    input:
    path filtFs_path
    path filtRs_path
    path filter_metadata

    // Define output channels
    output:
    path("dereps/*_F_*_dereped.RDS"), emit: dereps_F
    path("dereps/*_R_*_dereped.RDS"), emit: dereps_R

    script:
    """
    Rscript ${projectDir}/bin/dada2_dereplicate.R \
        --filtFs-path ${filtFs_path} \
        --filtRs-path ${filtRs_path} \
        --filter-metadata ${filter_metadata} \
        --dout dereps \
        --ncores ${task.cpus} \
        --log-level ${params.logLevel}
    """
}