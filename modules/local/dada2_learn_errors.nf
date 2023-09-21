// Define a Nextflow process for learning error models using DADA2
process DADA2_LEARN_ERRORS {

    label 'process_medium'

    // Define input channels
    input:
    path filt_1
    path filt_2
    path filter_metadata
    val maxConsist
    val randomize

    // Define output channels
    output:
    path("error_model/err_F_model.RDS"), optional: true, emit: error_model_F
    path("error_model/err_R_model.RDS"), optional: true, emit: error_model_R
    path("error_model/err_model.RDS"), optional: true, emit: error_model

    // Main script to run
    script:
    def filt_1_list = filt_1.join(' ')  
    def filt_2_list = filt_2.join(' ')  
    def randomize = randomize ? '--rand' : ''

    """
    Rscript ${projectDir}/bin/dada2_learn_errors.R \
        --filt-1 ${filt_1_list} \
        --filt-2 ${filt_2_list} \
        --ncores ${task.cpus} \
        --maxConsist ${maxConsist} \
        --dout error_model \
        ${randomize}
    """
}
