// workflows/quality_control.nf

include { PREPROCESS_COVERAGE } from '../modules/local/preprocess_coverage.nf'
include { POSTPROCESS_COVERAGE } from '../modules/local/postprocess_coverage.nf'
include { PUBLISH_COVERAGE } from '../modules/local/publish_coverage.nf'

workflow QUALITY_CONTROL {

    // Define inputs
    take:
    sample_coverage_files
    amplicon_coverage_files
    alleledata
    clusters

    main:

    // Assuming 'sampleFiles' and 'ampliconFiles' are your channels for individual files
    sample_combined = sample_coverage_files.collect()
    amplicon_combined = amplicon_coverage_files.collect()

    // Put all sample and amplicon coverage files into a single file with sample name column
    PREPROCESS_COVERAGE(
        sample_combined,
        amplicon_combined
    )
    
    // If postprocessing coverage is provided add dada2 and postprocessing coverage to the coverage files
    def postprocessing = (alleledata && clusters)
    
    if (postprocessing) {
        POSTPROCESS_COVERAGE(
            alleledata,
            clusters,
            PREPROCESS_COVERAGE.out.sample_coverage,
            PREPROCESS_COVERAGE.out.amplicon_coverage
        )
        sample_coverage_ch = POSTPROCESS_COVERAGE.out.postprocess_sample_coverage
        amplicon_coverage_ch = POSTPROCESS_COVERAGE.out.postprocess_amplicon_coverage
    } else {
        sample_coverage_ch   = PREPROCESS_COVERAGE.out.sample_coverage
        amplicon_coverage_ch = PREPROCESS_COVERAGE.out.amplicon_coverage
    }
    
    PUBLISH_COVERAGE(
        sample_coverage_ch,
        amplicon_coverage_ch
    )

    emit:
    sample_coverage = PUBLISH_COVERAGE.out.sample_coverage
    amplicon_coverage = PUBLISH_COVERAGE.out.amplicon_coverage
}
