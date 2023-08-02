/*
 * STEP - BUILD_RESISTANCE_TABLE
 * 
 */

process BUILD_RESISTANCE_TABLE {

  tag "$meta.id"
  label 'process_low'

  input:
  path alleledata
  path resmarkers
  path refseq

  output:
  path("alignments.pseudocigar.txt"), emit: pseudocigar

  script:
  """
  python3 ${projectDir}/bin/resistance_marker_module.py \
    --allele_data_path ${alleledata} \
    --res_markers_info_path ${resmarkers} \
    --refseq_path ${refseq} \
    --n-cores ${task.cpus}
  """
}