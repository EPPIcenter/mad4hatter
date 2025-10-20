/*
 * STEP - COLLAPSE_CONCATENATED_READS
 * Denoise the demultiplexed amplicon fastqs 
 */

process COLLAPSE_CONCATENATED_READS {

  label 'process_single'
  conda 'envs/postproc-env.yml'

  input:
  path clusters

  output:
  path 'clusters.concatenated.collapsed.txt', emit: clusters_concatenated_collapsed
  
  script:
  """
  Rscript ${projectDir}/bin/collapse_concatenated_reads.py \
    --clusters ${clusters} 
  """
}