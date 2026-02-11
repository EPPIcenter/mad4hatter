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

    def pool_name_mapping = [
        '1A': 'D1.1',
        '1B': 'R1.1',
        '2' : 'R2.1',
        '5' : 'R1.2', 
        'D1' : 'D1.1',
        'R1' : 'R1.2',
        'R2' : 'R2.1',
        'M1' : 'M1.1',
        'M2' : 'M2.1',
    ]

    def selectedPools_updated = []
    for (pool in selectedPools) {
        def paths = params.pool_options[pool.trim()]
        if (paths == null) {
            error "Pool '${pool}' not found in configuration."
        }
        amplicon_info_paths.add("$projectDir/$paths.amplicon_info_path")
        // if shorthand name replace pool with versioned name
        if (pool.trim() in pool_name_mapping) {
            selectedPools_updated.add(pool_name_mapping[pool.trim()])
        }
        else {
            selectedPools_updated.add(pool.trim())
        }
    }
    def amplicon_info_paths_str = amplicon_info_paths.join(' ')
    def selectedPools_str = selectedPools_updated.join(' ')

    BUILD_AMPLICON_INFO(selectedPools_str, amplicon_info_paths_str, "amplicon_info.tsv")

    // Emit the amplicon_info channel
    emit: amplicon_info_ch = BUILD_AMPLICON_INFO.out.amplicon_info
}
