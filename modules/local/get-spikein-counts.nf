// Process for getting spikein counts
process GET_SPIKEIN_COUNTS {
    tag "Getting Spikein Counts"
    label 'process_low'

    input:
    // The remaining reads to analyze; these fastqs should have
    // all panel reads removed at this point
    path (demultiplexed_spikeins_fastqs)

    output:
    path ('spikein_counts/*_final_spikein_counts.csv'), emit: spikein_counts
    path ('spikein_counts/*_multi_map_spikein_counts.csv'), emit: multimap

    publishDir(
        path: "${params.outDIR}/spikein_work",
        mode: 'copy',
        glob: '*.{txt,csv}'
    )

    script:
    """
    bash get-counts.sh \
        -c ${params.spikein_csv} \
        -d ${demultiplexed_spikeins_fastqs} \
        -o "spikein_counts" \
        -t ${task.cpus}
    """
}