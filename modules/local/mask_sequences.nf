// Dada2 Postprocessing
process MASK_SEQUENCES {

  conda 'pandoc'

  publishDir(
      path: "${params.outDIR}",
      mode: 'copy'
  )

  input:
  path masks
  path alignments
  val parallel
  val n_cores

  output:
  path "masked.alignments.txt", emit: masked_alignments

  script:
  def parallel = parallel ? '--parallel' : ''
  def n_cores = (parallel && n_cores > 0) ? "--n-cores ${n_cores}" : ''

  """
  Rscript ${projectDir}/bin/mask_sequences.R \
    --masks ${masks.join(' ')} \
    --alignments ${alignments} \
    ${parallel} \
    ${n_cores}

  """
}