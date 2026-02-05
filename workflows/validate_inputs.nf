nextflow.enable.dsl = 2

workflow VALIDATE_INPUTS {

    // Check if the workflow is valid
    def allowed_workflow_names = ['qc', 'complete', 'postprocessing']
    workflow_name = params.workflow_name?.toLowerCase()

    if (!allowed_workflow_names.contains(workflow_name)) {
        log.error "Invalid workflow specified: ${params.workflow_name}. Allowed workflows are: qc, complete, postprocessing."
        exit 1
    }

    // Check pools input
    check_pools()

    // Check params based on workflow
    if (workflow_name == 'complete' || workflow_name == 'qc') {
        check_readdir_presence()
        // Check sequencer input
        check_sequencer()
    } else if (workflow_name == 'postprocessing') {
        check_denoised_asvs_presence()
    }
    
    // Add warning if both params.refseq_fasta and params.genome are not null that params.refseq_fasta will be used
    if (params.refseq_fasta != null && params.genome != null) {
        log.warn "Both --refseq_fasta and --genome are specified. --refseq_fasta will be used and --genome will be ignored."
    }
}

// Helper function to check if pools parameter is provided
def check_pools() {
    if (params.pools == null) {
        log.error "`--pools` must be specified."
        exit 1
    }

    def legacy_pools = [
        '1A': 'D1 or D1.1',
        '1B': 'R1.1',
        '2' : 'R2 or R2.1',
        '5' : 'R1 or R1.2'
    ]
    def warnings = []

    def selectedPools = params.pools.toString().split(',')
    def invalidPools = []

    // Check all pools first
    for (pool in selectedPools) {
        def trimmedPool = pool.trim()
        def paths = params.pool_options[trimmedPool]
        if (paths == null) {
            invalidPools << trimmedPool
        } else {
            legacy_pools.each { legacy, current ->
                if (trimmedPool == legacy) {
                    warnings << "You have input a legacy pool name: ${legacy}. Current name would be ${current}."
                }
            }
        }
    }

    // Report all invalid pools at once
    if (!invalidPools.isEmpty()) {
        // Using pools outside of preset in configuration is okay as long as other required files are provided.
        def workflow_name = params.workflow_name?.toLowerCase()
        
        if (workflow_name == 'qc') {
            if (params.amplicon_info == null) {
                log.error "Pools were not found in configuration: ${invalidPools.join(', ')}. `--amplicon_info` must be specified when using bespoke pools."
                exit 1
            }
        } else if (workflow_name == 'complete' || workflow_name == 'postprocessing') {
            if (params.amplicon_info == null || (params.refseq_fasta == null && params.genome == null)) {
                log.error "Pools were not found in configuration: ${invalidPools.join(', ')}. `--amplicon_info` and either `--refseq_fasta` or `--genome` must be provided when using bespoke pools."
                exit 1
            }
        } else {
            log.warn "The following pools were not found in configuration: ${invalidPools.join(', ')}. Proceeding with bespoke pools."
        }
    }

    if (!warnings.isEmpty()) {
        warnings.each { warning -> log.warn warning }
    }
}

// Helper function to check if sequencer parameter is provided
def check_sequencer() {
    if (params.sequencer == null) {
        log.error "`--sequencer` must be specified."
        exit 1
    }
}

// Helper function to check if readDIR parameter is provided
def check_readdir_presence() {
    if (params.readDIR == null) {
        log.error "`--readDIR` must be specified but is missing."
        exit 1
    }
}

// Helper function to check if denoised_asvs parameter is provided
def check_denoised_asvs_presence() {
    if (params.denoised_asvs == null) {
        log.error "`--denoised_asvs` must be specified but is missing."
        exit 1
    }
}
