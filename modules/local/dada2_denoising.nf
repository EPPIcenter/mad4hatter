// Define a Nextflow module for DADA2 Denoising
process DADA2_DENOISING {

    label 'process_medium'

    // Define input channels
    input:
    path dereps_F
    path dereps_R
    path error_model_F
    path error_model_R
    val pool
    val band_size
    val omega_a
    val just_concatenate
    val use_quals
    val maxEE
    val self_consist
    val omega_c
    val verbose
    path ampliconFILE
    
    // Define output channels
    output:
    path("*seqtab.RDS"), emit: seqtab

    script:
    def concatenate = just_concatenate ? '--just-concatenate' : ''
    def verbose = verbose ? '--verbose' : ''
    def self_consist = self_consist ? '--self-consist' : ''

    """
    Rscript ${projectDir}/bin/dada2_denoising.R \
        --derep-1 ${dereps_F} \
        --derep-2 ${dereps_R} \
        --error-model-1 ${error_model_F} \
        --error-model-2 ${error_model_R} \
        --ncores ${task.cpus} \
        --pool ${pool} \
        --band-size ${band_size} \
        --omega-a ${omega_a} \
        --use-quals ${use_quals} \
        --maxEE ${maxEE} \
        --omega-c ${omega_c} \
        --ampliconFILE ${ampliconFILE} \
        --log-level ${params.logLevel} \
        ${verbose} \
        ${concatenate} \
        ${self_consist}
    """
}