/*
 * STEP - BUILD_PSEUDOCIGAR
 * 
 */

process BUILD_PSEUDOCIGAR {

  label 'process_medium'
  conda 'envs/postproc-env.yml'

  input:
  path alignments
  val output_suffix

  output:
  path("alignments.pseudocigar${output_suffix}.txt"), emit: pseudocigar

  script:
  """
  Rscript ${projectDir}/bin/build_pseudocigar.R \
    --alignments ${alignments} \
    --output_suffix ${output_suffix} \
    --ncores ${task.cpus}
  """
}