nextflow.enable.dsl=2

include { BUILD_AMPLICON_INFO } from '../modules/local/build_resources.nf'


workflow GENERATE_AMPLICON_INFO {
    // Check if pools parameter is provided
    if (params.pools == null) {
        error "Please specify the pools using --pools"
    }
    // Split the pools parameter into a list
    def selectedPools = params.pools.toString().split(',')

    // List resource paths 
    def amplicon_info_paths = []

    for (pool in selectedPools) {
        def paths = params.pool_options[pool.trim()]
        if (paths == null) {
            error "Pool '${pool}' not found in configuration."
        }
        amplicon_info_paths.add("$projectDir/$paths.amplicon_info_path")
    }
    def amplicon_info_paths_str = amplicon_info_paths.join(' ')
    def selectedPools_str = selectedPools.join(' ')

    BUILD_AMPLICON_INFO(selectedPools_str, amplicon_info_paths_str, "amplicon_info.tsv")

    // Emit the amplicon_info channel
    emit: amplicon_info_ch = BUILD_AMPLICON_INFO.out.amplicon_info
}
