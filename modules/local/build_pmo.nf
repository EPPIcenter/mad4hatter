/*
 * STEP - BUILD_PMO
 * Combine the pipeline outputs into a Portable Microhaplotype Object (PMO) and validate it
 */

process BUILD_PMO {

  label 'process_single'
  conda 'envs/pmo-env.yml'

  publishDir(path: "${params.outDIR}", pattern: "pmo.json.gz", mode: 'copy')

  input:
  path allele_data
  path amplicon_info
  val panel_name
  path specimen_key

  output:
  path("pmo.json.gz"), emit: pmo

  script:
  def specimen_key_arg = specimen_key ? "--specimen_key ${specimen_key}" : ''
  """
  \${PMO_PYTHON:-python3} ${projectDir}/bin/build_pmo.py \
    --allele_data ${allele_data} \
    --amplicon_info ${amplicon_info} \
    --panel_name ${panel_name} \
    --bioinformatics_run_name ${workflow.runName} \
    --output pmo.json.gz \
    ${specimen_key_arg}
  """
}
