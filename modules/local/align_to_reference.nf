/*
 * STEP - ALIGN-TO-REFERENCE
 * 
 */

process ALIGN_TO_REFERENCE {

    input:
    path clusters
    path refseq_fasta
    path amplicon_info
    val parallel
    val n_cores

    output:
    path("alignments.txt"), emit: alignments

    script:
    def parallel = parallel ? '--parallel' : ''
    def n_cores = (parallel && n_cores > 0) ? "--n-cores ${n_cores}" : ''

    """
    Rscript ${projectDir}/bin/align_to_reference.R \
      --clusters ${clusters} \
      --refseq-fasta ${refseq_fasta} \
      --amplicon-table ${amplicon_info} \
      ${parallel} \
      ${n_cores}
    """
}