// workflows/quality_control.nf

include { PREPROCESS_COVERAGE } from '../modules/local/preprocess_coverage.nf'
include { POSTPROCESS_COVERAGE } from '../modules/local/postprocess_coverage.nf'
include { QUALITY_REPORT } from '../modules/local/quality_report.nf'

workflow QUALITY_CONTROL {

    // Define inputs
    take:
    amplicon_info
    sample_coverage_files
    amplicon_coverage_files
    alleledata
    clusters

    main:

    // Assuming 'sampleFiles' and 'ampliconFiles' are your channels for individual files
    sample_combined = sample_coverage_files.collect()
    amplicon_combined = amplicon_coverage_files.collect()

    // Initial Preprocessing
    PREPROCESS_COVERAGE(
        sample_combined,
        amplicon_combined
    )

    // If postprocessing coverage is provided, run the postprocessing workflow
    def postprocessing = (alleledata && clusters)
    sample_coverage_ch = postprocessing ? 
        POSTPROCESS_COVERAGE(
            alleledata,
            clusters,
            PREPROCESS_COVERAGE.out.sample_coverage,
            PREPROCESS_COVERAGE.out.amplicon_coverage
        ).postprocess_sample_coverage :
        PREPROCESS_COVERAGE.out.sample_coverage

    amplicon_coverage_ch = postprocessing ? 
        POSTPROCESS_COVERAGE.out.postprocess_amplicon_coverage : 
        PREPROCESS_COVERAGE.out.amplicon_coverage

    // Reporting
    QUALITY_REPORT(
        sample_coverage_ch,
        amplicon_coverage_ch,
        amplicon_info
    )
}
