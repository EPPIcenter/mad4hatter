/*
 * STEP - BUILD_RESISTANCE_TABLE
 * 
 */

process BUILD_RESISTANCE_TABLE {

  label 'process_medium'
  conda 'envs/resmarker-env.yml'

  publishDir(
    path: "${params.outDIR}/resistance_marker_module",
    mode: 'copy'
  )

  input:
  path alleledata
  path resmarkers
  path refseq

  output:
  path ('resmarker_table.txt'), emit: resmarkers
  path ('resmarker_microhap_table.txt'), emit: microhaps
  path ('resmarker_new_mutations.txt'), emit: new_mutations

  script:
  """
  python3 ${projectDir}/bin/resistance_marker_module.py \
    --allele_data_path ${alleledata} \
    --res_markers_info_path ${resmarkers} \
    --refseq_path ${refseq} \
    --n-cores ${task.cpus} 
  """
}