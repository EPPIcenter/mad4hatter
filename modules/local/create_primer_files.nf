/*
 * STEP - CREATE_PRIMER_FILES
 * Prepare the primer files from the given amplicon_info file
 */

process CREATE_PRIMER_FILES {

    tag "$meta.id"
    label 'process_single'

    input:
    path amplicon_info

    output:
    path("fwd_primers.fasta"), emit: fwd_primers
    path("rev_primers.fasta"), emit: rev_primers

    shell:
    """
    bash create_primer_files.sh -a ${amplicon_info} -f fwd_primers.fasta -r rev_primers.fasta
    """
}