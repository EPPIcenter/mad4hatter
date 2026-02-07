/*
 * STEP - CUTADAPT
 * Prepare the primer files from the given amplicon_info file
 */


process CUTADAPT {

  tag "$pair_id"
  label 'process_low'
  conda 'envs/cutadapt-env.yml'

    input:
    file fwd_primers
    file rev_primers
    tuple val(pair_id), file(reads)
    val cutadapt_minlen
    val gtrim
    val allowed_errors

    output:
    path("*.SAMPLEsummary.txt"), emit: sample_summary
    path("*.AMPLICONsummary.txt"), emit: amplicon_summary
    path('demultiplexed_fastqs'), emit: demultiplexed_fastqs

    script:
    """
    bash cutadapt_process.sh \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -r ${rev_primers} \
        -f ${fwd_primers} \
        -m ${cutadapt_minlen} \
        -g ${gtrim} \
        -e ${allowed_errors} \
        -c ${task.cpus} \
        -o demultiplexed_fastqs
    """
}