/*
 * STEP - BUILD_PSEUDOCIGAR
 * 
 */

process BUILD_PSEUDOCIGAR {

  label 'process_medium'
  conda 'envs/postproc-env.yml'

  input:
  path alignments
  val output_name 

  output:
  path("alignments.pseudocigar*.txt"), emit: pseudocigar

  script:
  """
  Rscript ${projectDir}/bin/build_pseudocigar.R \
    --alignments ${alignments} \
    --output_suffix ${output_name} \
    --ncores ${task.cpus}
  """
}