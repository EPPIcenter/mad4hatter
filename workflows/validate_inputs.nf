nextflow.enable.dsl = 2

workflow VALIDATE_INPUTS {

    // Check if the workflow is valid
    def allowed_workflows = ['qc', 'complete', 'postprocessing']
    workflow = params.workflow?.toLowerCase()

    if (!allowed_workflows.contains(workflow)) {
        log.error "Invalid workflow specified: ${params.workflow}. Allowed workflows are: qc, complete, postprocessing."
        exit 1
    }

    // Check pools input
    check_pools()

    // Check sequencer input
    check_sequencer()

    // Check params based on workflow
    if (workflow == 'complete' || workflow == 'qc') {
        check_readdir_presence()
    } else if (workflow == 'postprocessing') {
        check_denoised_asvs_presence()
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
