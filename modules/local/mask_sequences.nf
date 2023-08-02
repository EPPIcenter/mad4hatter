// Dada2 Postprocessing
process MASK_SEQUENCES {

  tag "$meta.id"
  label 'process_low'

  conda 'pandoc'

  input:
  path masks
  path alignments

  output:
  path "masked.alignments.txt", emit: masked_alignments

  script:
  def n_cores = "${task.cpus}" ? "--n-cores ${task.cpus}" : ''

  """
  Rscript ${projectDir}/bin/mask_sequences.R \
    --masks ${masks.join(' ')} \
    --alignments ${alignments}

  """
}