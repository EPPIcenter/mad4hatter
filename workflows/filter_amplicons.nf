workflow FILTER_AMPLICONS {

    take:
    demultiplexed_fastqs

    main:

    // Step 1: dada2_preprocessing
    DADA2_FILTERING(
        demultiplexed_fastqs,
        params.minLen,
        params.maxN,
        params.rm_phix,
        params.compress,
        params.matchIDs,
        params.maxEE_R1,
        params.truncQ_R1,
        params.trimRight_R1,
        params.trimLeft_R1,
        params.maxEE_R2,
        params.truncQ_R2,
        params.trimRight_R2,
        params.trimLeft_R2
    )

    emit:
    dada_filtFs_ch = DADA2_FILTERING.out.filtFs
    dada_filtRs_ch = DADA2_FILTERING.out.filtRs
    dada_filter_metadata_ch = DADA2_FILTERING.out.filter_metadata
}