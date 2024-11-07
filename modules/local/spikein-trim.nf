process SPIKEIN_TRIM {

    label 'process_low'

    input:
    path fwd_primers
    path rev_primers
    path unknown_fastqs

    output:
    path("spikein_trim/demultiplexed/*"), emit: spikeins_demultiplexed

    shell:
    """
    bash spikein_trim.sh \
        -f ${unknown_fastqs} \
        -1 ${fwd_primers} \
        -2 ${rev_primers} \
        -r "spikein_trim" \
        -t ${task.cpus}
    """
}