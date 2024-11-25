// Process for generating plots
process PLOT_SPIKEIN_COUNTS {
    tag "Generating Plots"
    label 'process_single'

    publishDir(
        path: "${params.outDIR}/quality_report",
        mode: 'copy',
        pattern: 'contamination_report.pdf'
    )

    publishDir(
        path: "${params.outDIR}/quality_report",
        mode: 'copy',
        pattern: 'spikein_amplicon_ratio.csv'
    )

    input:
    path (counts_files, stageAs: "counts_files?")
    path (expected_spikein)
    path (spikein_info)
    path (amplicon_coverage)

    output:
    path 'contamination_report.pdf', emit: contamination_report
    path 'spikein_amplicon_ratio.csv', emit: spikein_amplicon_ratio

    script:
    """
    Rscript ${projectDir}/bin/spikein_detection_heatmap.R \
        --input ${counts_files} \
        --expected ${expected_spikein} \
        --spikein-info ${spikein_info} \
        --output contamination_report.pdf \
        --amplicon-coverage ${amplicon_coverage}
    """
}