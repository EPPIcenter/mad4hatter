/*
 * STEP - DADA2_ANALYSIS
 * Denoise the demultiplexed amplicon fastqs 
 */

process DADA2_ANALYSIS {

  label 'process_high'
  conda 'envs/dada2-env.yml'

  publishDir(
    path: "${params.outDIR}/raw_dada2_output",
    mode: 'copy',
    pattern: '{dada2.clusters.txt,*.RDS}'
  )

  input:
  path (demultiplexed_fastqs, stageAs: "demultiplexed_fastqs?")
  path amplicon_info
  val pool
  val band_size
  val omega_a
  val maxEE
  val just_concatenate

  output:
  path 'dada2.clusters.txt', emit: dada2_clusters
  path '*.RDS', emit: rds_files  // Make sure to emit this output

  
  script:
  
  def concatenate = just_concatenate ? '--concat-non-overlaps' : ''

  """
  Rscript ${projectDir}/bin/dada_overlaps.R \
    --trimmed-path ${demultiplexed_fastqs} \
    --ampliconFILE ${amplicon_info} \
    --pool ${params.pool} \
    --band-size ${params.band_size} \
    --omega-a ${params.omega_a} \
    --maxEE ${params.maxEE} \
    --cores ${task.cpus} \
    ${concatenate} 
  """
}