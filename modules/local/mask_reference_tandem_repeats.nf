process MASK_REFERENCE_TANDEM_REPEATS {

  tag "$meta.id"
  label 'process_single'

  conda 'envs/trf-env.yml'

  input:
  path refseq_fasta
  val min_score
  val max_period

  output:
  path "*.mask", emit: masked_fasta

  script:
  """
  trf ${refseq_fasta} 2 7 7 80 10 ${min_score} ${max_period} -h -m
  """
}