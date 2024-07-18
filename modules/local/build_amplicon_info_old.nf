// // /*
// //  * STEP - BUILD_AMPLICON_INFO
// //  * 
// //  */


// // params.targets = ['1A','1B']
// // params.amplicon_info_paths = ['resources/1A/1A_amplicon_info.tsv','resources/1B/1B_amplicon_info.tsv']

// // process addTarget {

// //     input:
// //     file amplicon_info_path from amplicon_info_paths
// //     val target

// //     output:
// //     file("${amplicon_info_path.baseName}_with_target.txt") into concat_amplicon_info

// //     script:
// //     """
// //     awk -v target=${target} '{print \$0, "\t", target}' ${amplicon_info_path} > ${amplicon_info_path.baseName}_with_target.txt
// //     """
// // }

// // process concatenateAmpliconInfo {

// //     input:
// //     file concat_amplicon_info_ch

// //     output:
// //     file('concatenated_amplicon_info.txt')

// //     script:
// //     """
// //     cat ${concat_amplicon_info_ch.collect()} | sort -u > concatenated_amplicon_info.txt
// //     """
// // }

// // workflow {
// //     // Convert single values to lists if necessary
// //     targets = params.targets instanceof List ? params.targets : [params.targets]
// //     amplicon_info_paths = params.amplicon_info_paths instanceof List ? params.amplicon_info_paths : [params.amplicon_info_paths]
    
// //     // Combine into one channel 
// //     combinedChannel = Channel.from(targets).merge(Channel.from(amplicon_info_paths))
// //     combinedChannel.view()

// //     // Add target column to each amplicon_info file
// //     addTarget(input: combinedChannel)
// //     // Concatenate all the files and remove duplicates
// //     // concatenateAmpliconInfo(concat_amplicon_info_ch)
// // }
// // nextflow.enable.dsl=2

// // process annotate_csv {
// //     input:
// //     tuple val(pool), path(csv_file)
    
// //     output:
// //     path "${pool}_annotated.csv"
    
// //     script:
// //     """
// //     python process_amplicons.py --pool ${pool} --input ${csv_file} --output ${pool}_annotated.csv
// //     """
// // }

// // process concatenate_csv {
// //     input:
// //     path annotated_csvs
    
// //     output:
// //     path "concatenated_amplicons.csv"
    
// //     script:
// //     """
// //     python process_amplicons.py --concatenate --inputs ${annotated_csvs} --output concatenated_amplicons.csv
// //     """
// // }

// // params.pools = ["1A", "1B", "2"]
// // params.csv_files = ["/Users/kmurie/Documents/git_projects/mad4hatter-dev/mad4hatter/resources/1A/1A_amplicon_info.tsv", "/Users/kmurie/Documents/git_projects/mad4hatter-dev/mad4hatter/resources/1B/1B_amplicon_info.tsv", "/Users/kmurie/Documents/git_projects/mad4hatter-dev/mad4hatter/resources/2/2_amplicon_info.tsv"]
    

// // workflow {
// //     pools = params.pools instanceof List ? params.pools : [params.pools]

// //     // 
// //     combinedChannel = Channel.from(params.pools).merge(Channel.from(params.csv_files))
// //     combinedChannel.view()

// //     annotate_csv(combinedChannel)
// //     // concatenate_csv(annotated_csvs.flatten())
// // }

// // Pipeline takes either list or single value for targets. 
// // For pool read paths to amplicon info from config 
// // Generate compiled amplicon info
// // For pool read paths to resmarker from config 
// // Generate compiled resmarker table 


// nextflow.enable.dsl=2

// process annotate_csv {
//     input:
//     tuple val(pool), path(csv_file)

//     output:
//     path "${pool}_annotated.csv"

//     script:
//     """
//     python process_amplicons.py --pool ${pool} --input ${csv_file} --output ${pool}_annotated.csv
//     """
// }

// process concatenate_csv {
//     input:
//     path annotated_csvs

//     output:
//     path "concatenated_amplicons.csv"

//     script:
//     """
//     python process_amplicons.py --concatenate --inputs ${annotated_csvs.join(' ')} --output concatenated_amplicons.csv
//     """
// }

// workflow {
//     // Get the list of pools from the command-line argument
//     def selected_pools = params.pools.split(',')
//     print(selected_pools)
//     // // Load the configuration
//     // def pool_configs = [:]
//     // for (pool in selected_pools) {
//     //     def pool_conf = config.pools[pool.trim()]
//     //     if (pool_conf) {
//     //         pool_configs[pool.trim()] = pool_conf
//     //     } else {
//     //         error "Pool '${pool}' not found in the configuration"
//     //     }
//     // }

//     // // Create a list of pool names and their corresponding amplicon_info_path
//     // def amplicon_info_paths = pool_configs.collect { pool, conf ->
//     //     [pool, file(conf.amplicon_info_path)]
//     // }
//     // amplicon_info_paths.value()
//     // // Create a channel with tuples of pool names and corresponding CSV files
//     // annotate_results = Channel
//     //     .from(amplicon_info_paths)
//     //     .map { pool, csv_file -> tuple(pool, amplicon_info_path) }
//     //     .set { annotated_csvs }

//     // // Run the annotation process
//     // annotate_csv(annotated_csvs)

//     // // Run the concatenation process
//     // concatenate_csv(annotated_csvs.flatten())
// }


// Function to get pool paths from the config
def get_pool_paths(poolName) {
    def poolConfig = params.pools_config[poolName]
    return [
        amplicon_info_path: poolConfig.amplicon_info_path,
        resistance_marker_table_path: poolConfig.resistance_marker_table_path,
        refseq_path: poolConfig.refseq_path
    ]
}

workflow {
    // Check if pools parameter is provided
    if (params.pools == null) {
        error "Please specify the pools using --pools"
    }

    // // Split the pools parameter into a list
    // def selectedPools = params.pools.split(',')

    // // Initialize lists to hold the paths
    // def amplicon_info_paths = []
    // def resistance_marker_table_paths = []
    // def refseq_paths = []

    // // Populate the lists with paths from the config
    // for (pool in selectedPools) {
    //     def paths = get_pool_paths(pool)
    //     amplicon_info_paths << paths.amplicon_info_path
    //     resistance_marker_table_paths << paths.resistance_marker_table_path
    //     refseq_paths << paths.refseq_path
    // }

    // // Print the lists to verify
    // println "Amplicon Info Paths: ${amplicon_info_paths}"
    // println "Resistance Marker Table Paths: ${resistance_marker_table_paths}"
    // println "RefSeq Paths: ${refseq_paths}"

    // // Now you can use these lists in your workflow
    // // For example, passing them to a process
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
