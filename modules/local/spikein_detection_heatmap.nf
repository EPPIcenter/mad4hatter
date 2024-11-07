// Process for generating plots
process SPIKEIN_DETECTION {
    tag "Generating Plots"
    publishDir path: results_path(''), mode: 'copy'

    input:
    path counts_file from counts

    output:
    path results_path('spikein_counts_heatmap.png')

    script:
    """
    Rscript spikein_detection_heatmap.R \
        --input ${counts_file} \
        --output ${}
    """
}

