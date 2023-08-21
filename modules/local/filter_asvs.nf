/*
 * STEP - ASV FITLERING
 * Prepare the primer files from the given amplicon_info file
 */


process FILTER_ASVS {

    label 'process_low'

    input:
    path alignments

    output:
    path("filtered.alignments.txt"), emit: filtered_alignments_ch

    script:
    """
    bash filter_asv_process.sh \
        -i ${alignments} \
        -o filtered.alignments.txt \
        -t ${params.alignment_threshold}
    """
}