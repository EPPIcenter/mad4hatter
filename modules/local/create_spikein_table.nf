/*
 * STEP - ALIGN-TO-REFERENCE
 * 
 */

process CREATE_SPIKEIN_TABLE {

  label 'process_medium'

  input:
  path spikein_fastas
  
  output:
  path("spikeins_data.txt"), emit: spikein_data

  publishDir(
    path: "${params.outDIR}/stats",
    mode: 'copy',
    pattern: '*.txt'
  )

  script:
  """
  Rscript ${projectDir}/bin/create_spikein_table.R \
    --trimmed-path ${spikein_fastas} \
    --maxEE ${params.maxEE}
  """
}