// DADA2_PREPROCESSING Process

// Filters and dereplicates reads using DADA2

process DADA2_PREPROCESSING {
    // Define CPU and memory requirements
    label 'process_low'

    // Define input channels
    input:
    path(trimmed_path)

    val(minLen)
    val(maxN)
    val(rm_phix)
    val(compress)
    val(matchIDs)
    val(maxEE_R1)
    val(truncQ_R1)
    val(trimRight_R1)
    val(trimLeft_R1)
    val(maxEE_R2)
    val(truncQ_R2)
    val(trimRight_R2)
    val(trimLeft_R2)

    // Define output channels
    output:
    path("filtered_${trimmed_path.baseName}/*_F_filt.fastq.gz"), emit: filtFs
    path("filtered_${trimmed_path.baseName}/*_R_filt.fastq.gz"), emit: filtRs
    path("filtered_${trimmed_path.baseName}/*_filter_metadata.RDS"), emit: filter_metadata

    // Main script to run
    script:
    def verbose = params.verbose ? "--verbose" : ""
    
    """
    Rscript ${projectDir}/bin/dada2_preprocessing.R \
      --trimmed-path ${trimmed_path} \
      --maxN ${maxN} \
      --rm-phix ${rm_phix} \
      --compress ${compress} \
      --ncores ${task.cpus} \
      --minLen ${minLen} \
      --maxEE_R1 ${maxEE_R1} \
      --truncQ_R1 ${truncQ_R1} \
      --trimRight_R1 ${trimRight_R1} \
      --trimLeft_R1 ${trimLeft_R1} \
      --maxEE_R2 ${maxEE_R2} \
      --truncQ_R2 ${truncQ_R2} \
      --trimRight_R2 ${trimRight_R2} \
      --trimLeft_R2 ${trimLeft_R2} \
      --matchIDs ${matchIDs} \
      ${verbose}
    """
}
