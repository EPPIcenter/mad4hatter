// Process for generating plots
process PLOT_SPIKEIN_COUNTS {
    tag "Generating Plots"
    label 'process_single'

    publishDir(
        path: "${params.outDIR}/quality_report",
        mode: 'copy',
        pattern: 'contamination_report.pdf'
    )

    input:
    path (counts_files, stageAs: "counts_files?")
    path (expected_spikein)
    path (spikein_info)
    path (amplicon_coverage)

    output:
    path 'contamination_report.pdf', emit: contamination_report

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