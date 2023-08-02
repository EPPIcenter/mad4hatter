/*
 * STEP - ALIGN-TO-REFERENCE
 * 
 */

process ALIGN_TO_REFERENCE {
  
  tag "$meta.id"
  label 'process_low'

  input:
  path clusters
  path refseq_fasta
  path amplicon_info
  
  output:
  path("alignments.txt"), emit: alignments

  script:
  def n_cores = "${task.cpus}" ? "--n-cores ${task.cpus}" : ''

  """
  Rscript ${projectDir}/bin/align_to_reference.R \
    --clusters ${clusters} \
    --refseq-fasta ${refseq_fasta} \
    --amplicon-table ${amplicon_info} \
    ${n_cores}
  """
}