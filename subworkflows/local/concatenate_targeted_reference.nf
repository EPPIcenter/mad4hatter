nextflow.enable.dsl=2

include { BUILD_TARGETED_REFERENCE } from '../../modules/local/build_resources.nf'

workflow CONCATENATE_TARGETED_REFERENCE {
    // Check if pools parameter is provided
    if (params.pools == null) {
        error "Please specify the pools using --pools"
    }
    // Split the pools parameter into a list
    def selectedPools = params.pools.toString().split(',')

    // List resource paths 
    def targeted_reference_paths = []

    for (pool in selectedPools) {
        def paths = params.pool_options[pool.trim()]
        if (paths == null) {
            error "Pool '${pool}' not found in configuration."
        }
        targeted_reference_paths.add("$projectDir/$paths.targeted_reference_path")
    }
    def targeted_reference_paths_str = targeted_reference_paths.join(' ')
    BUILD_TARGETED_REFERENCE(targeted_reference_paths_str, "reference.fasta")

    // Emit the amplicon_info channel
    emit: reference_fasta = BUILD_TARGETED_REFERENCE.out.reference_fasta
}