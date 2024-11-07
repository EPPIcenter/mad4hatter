// Process for generating plots
process SPIKEIN_DETECTION {
    tag "Generating Plots"
    publishDir(
        path: "${params.outDIR}/quality_report",
        mode: 'copy'
    )

    input:
    path (counts_file)

    output:
    path ('*.pdf')

    script:
    """
    Rscript spikein_detection_heatmap.R \
        --input ${counts_file} \
        --expected ${params.expected_spikein} \
        --spikein-info ${params.spikein_info} \
        --output contamination_report.pdf
    """
}

