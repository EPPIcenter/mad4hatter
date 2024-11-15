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

    // Check params based on workflow
    if (workflow == 'complete' || workflow == 'qc') {
        check_readdir_presence()
        // Check sequencer input
        check_sequencer()
        // Check spike-in metadata
        check_spikein_metadata()

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

    def legacy_pools = [
        '1A': 'D1 or D1.1',
        '1B': 'R1.1',
        '2' : 'R2 or R2.1',
        '5' : 'R1 or R1.2'
    ]
    def warnings = []

    legacy_pools.each { legacy, current ->
        if (params.pools.toString().contains(legacy)) {
            warnings << "You have input a legacy pool name: ${legacy}. Current name would be ${current}."
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

// Helper function to check if spike-in metadata is
// available and complete.
def check_spikein_metadata() {

    // Map parameter names to their values
    def spikein_metadata_params = [
        "expected_spikein": params.expected_spikein,
        "spikein_info": params.spikein_info,
        "spikein_csv": params.spikein_csv,
        "spikein_primers": params.spikein_primers
    ]
    
    // Count null values
    def null_values = spikein_metadata_params.values().count { it == null }

    // If some files are present but others are missing, flag this as an error
    if (null_values != 4 && null_values != 0) {
        def missing_files = spikein_metadata_params.findAll { it.value == null }.keySet()
        log.error "Some spike-in metadata files are missing. Missing files are: ${missing_files.join(',')}"
        exit 1
    }
}

