// Process for generating plots
process PLOT_SPIKEIN_COUNTS {
    tag "Generating Plots"
    label 'process_single'

    publishDir(
        path: "${params.outDIR}/plots",
        mode: 'copy',
        pattern: '*.png'
    )

    publishDir(
        path: "${params.outDIR}/stats",
        mode: 'copy',
        pattern: '*.csv'
    )

    input:
    path (counts_files, stageAs: "counts_files?")

    output:
    path 'spikein_counts_heatmap.png', emit: spikein_counts_heatmap
    path 'spikein_counts_data.csv', emit: spikein_counts_data

    script:
    """
    Rscript ${projectDir}/bin/spikein_detection_heatmap.R \
        --input ${counts_files} \
        --output spikein_counts_heatmap.png
    """
}