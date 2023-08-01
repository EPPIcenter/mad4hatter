/*
 * STEP - COLLAPSE_CONCATENATED_READS
 * Denoise the demultiplexed amplicon fastqs 
 */

process COLLAPSE_CONCATENATED_READS {

  input:
  path clusters

  output:
  path 'clusters.concatenated.collapsed.txt', emit: clusters_concatenated_collapsed
  
  script:
  """
  Rscript ${projectDir}/bin/collapse_concatenated_reads.R \
    --clusters ${clusters} 
  """
}