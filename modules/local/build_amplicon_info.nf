// build_amplicon_info.nf
workflow {
    // Check if pools parameter is provided
    if (params.pools == null) {
        error "Please specify the pools using --pools"
    }

    // Split the pools parameter into a list
    def selectedPools = params.pools.split(',')
    print(pool_options)
    // Initialize lists to hold the paths
    def amplicon_info_paths = []
    def resistance_marker_table_paths = []
    def refseq_paths = []

    // Populate the lists with paths from the config
    for (pool in selectedPools) {
        def paths = params.pool_options[pool.trim()]
        if (paths == null) {
            error "Pool '${pool}' not found in configuration."
        }
        amplicon_info_paths << paths.amplicon_info_path
        resistance_marker_table_paths << paths.resistance_marker_table_path
        refseq_paths << paths.refseq_path
    }

    // Print the lists to verify
    println "Amplicon Info Paths: ${amplicon_info_paths}"
    println "Resistance Marker Table Paths: ${resistance_marker_table_paths}"
    println "RefSeq Paths: ${refseq_paths}"

    // // Example process using these lists
    // process buildAmpliconInfo {
    //     input:
    //     val amplicon_info_paths
    //     val resistance_marker_table_paths
    //     val refseq_paths

    //     script:
    //     """
    //     echo "Amplicon Info Paths: ${amplicon_info_paths.join(', ')}"
    //     echo "Resistance Marker Table Paths: ${resistance_marker_table_paths.join(', ')}"
    //     echo "RefSeq Paths: ${refseq_paths.join(', ')}"
    //     """
    // }
}
