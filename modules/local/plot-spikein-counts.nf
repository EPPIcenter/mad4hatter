// Process for generating plots
process PLOT_SPIKEIN_COUNTS {
    tag "Generating Plots"
    label 'process_single'

    publishDir(
        path: "${params.outDIR}/quality_report",
        mode: 'copy',
        pattern: '*.pdf'
    )

    input:
    path (counts_files, stageAs: "counts_files?")

    output:
    path 'contamination_report.pdf', emit: contamination_report

    script:
    """
    Rscript ${projectDir}/bin/spikein_detection_heatmap.R \
        --input ${counts_files} \
        --expected ${params.expected_spikein} \
        --spikein-info ${params.spikein_info} \
        --output contamination_report.pdf
    """
}