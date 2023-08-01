/*
 * STEP - CUTADAPT
 * Prepare the primer files from the given amplicon_info file
 */


process CUTADAPT {

    input:
    file fwd_primers
    file rev_primers
    tuple val(pair_id), file(reads)
    val cutadapt_minlen
    val sequencer
    val allowed_errors
    val cores

    output:
    path("*.SAMPLEsummary.txt"), emit: sample_summary
    path("*.AMPLICONsummary.txt"), emit: amplicon_summary
    path('demuliplexed_fastqs'), emit: demultiplexed_fastqs

    script:
    """
    bash cutadapt_process.sh \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -r ${rev_primers} \
        -f ${fwd_primers} \
        -m ${cutadapt_minlen} \
        -s ${sequencer} \
        -e ${allowed_errors} \
        -c ${cores} \
        -o demuliplexed_fastqs
    """
}