// Process for generating plots
process PLOT_SPIKEIN_COUNTS {
    tag "Generating Plots"
    label 'process_single'

    publishDir(
        path: "${params.outDIR}/contamination_report",
        mode: 'copy',
        pattern: 'contamination_report.pdf'
    )

    input:
    path (counts_files, stageAs: "counts_files?")
    path (alleledata)

    output:
    path 'contamination_report.pdf', emit: contamination_report

    script:
    // Rf. https://nextflow-io.github.io/patterns/optional-input/
    def alleledata_table = alleledata.name != "NO_FILE" ? "--alleledata ${alleledata}" : ""

    """
    Rscript ${projectDir}/bin/spikein_detection_heatmap.R \
        --input ${counts_files} \
        --expected ${params.expected_spikein} \
        --spikein-info ${params.spikein_info} \
        --output contamination_report.pdf \
        ${alleledata_table}
    """
}