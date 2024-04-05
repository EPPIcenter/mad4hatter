/*
 * STEP - ALIGN-TO-REFERENCE
 * 
 */

process ALIGN_TO_REFERENCE {

  label 'process_high'
  conda 'envs/postproc-env.yml'

  input:
  path clusters
  path refseq_fasta
  path amplicon_info
  
  output:
  path("alignments.txt"), emit: alignments

  script:
  
  """
  Rscript ${projectDir}/bin/align_to_reference.R \
    --clusters ${clusters} \
    --refseq-fasta ${refseq_fasta} \
    --amplicon-table ${amplicon_info} \
    --n-cores ${task.cpus}
  """
}