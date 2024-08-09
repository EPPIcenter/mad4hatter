nextflow.enable.dsl = 2

workflow VALIDATE_INPUTS {

    // Check if the workflow is valid
    def allowed_workflows = ['qc', 'complete', 'postprocessing']
    params.workflow = params.workflow?.toLowerCase()

    if (!allowed_workflows.contains(params.workflow)) {
        log.error "Invalid workflow specified: ${params.workflow}. Allowed workflows are: qc, complete, postprocessing."
        exit 1
    }

    // Check pools input
    check_pools()

    // Check sequencer input
    check_sequencer()

    // Check params based on workflow
    if (params.workflow == 'complete' || params.workflow == 'qc') {
        check_readdir_presence()
    } else if (params.workflow == 'postprocessing' && params.denoised_asvs != null) {
        check_denoised_asvs_presence()
    }

    // Define valid parameters and their valid values
    def valid_params = [
        'pools', 'workflow', 'sequencer', 'readDIR',
        'denoised_asvs', 'outDIR', 'profile', 'omega_a', 'dada2_pool', 'maxEE',
        'concat_non_overlaps', 'refseq_fasta', 'genome',
        'homopolymer_threshold', 'trf_min_score', 'trf_max_period', 'config'
    ]
    
    def valid_param_with_options = [
        'workflow'   : ['qc', 'complete', 'postprocessing'],
        'sequencer'  : ['miseq', 'nextseq'],
        'profile'    : ['sge', 'apptainer', 'docker', 'conda', 'slurm', 'mamba'],
        'dada2_pool' : ['pseudo', 'true', 'false']
    ]

    // Check if user-provided parameters are valid
    params.each { k, v ->
        // Check if the parameter is recognized
        if (!valid_params.contains(k)) {
            log.error "Unrecognized parameter: ${k}"
            exit 1
        }

        // Check if the parameter value is valid (if there are specific options)
        if (valid_param_with_options.containsKey(k) && !valid_param_with_options[k].contains(v)) {
            log.error "Invalid value for parameter: ${k}. Allowed values are: ${valid_param_with_options[k].join(', ')}"
            exit 1
        }
    }
}

// Helper function to check if pools parameter is provided
def check_pools() {
    if (params.pools == null) {
        log.error "`--pools` must be specified."
        exit 1
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
