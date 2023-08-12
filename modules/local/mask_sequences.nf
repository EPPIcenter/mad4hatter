// Dada2 Postprocessing
process MASK_SEQUENCES {

  label 'process_medium'

  conda 'pandoc'

  input:
  path masks
  path alignments

  output:
  path "masked.alignments.txt", emit: masked_alignments

  script:
  
  """
  Rscript ${projectDir}/bin/mask_sequences.R \
    --masks ${masks.join(' ')} \
    --alignments ${alignments} \
    --n-cores ${task.cpus}

  """
}