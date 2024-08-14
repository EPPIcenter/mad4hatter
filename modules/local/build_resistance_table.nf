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
  path alignment_data
  path resmarkers
  path refseq

  output:
  path ('resmarker_table.txt'), emit: resmarkers
  path ('resmarker_table_by_locus.txt'), emit: resmarkers_by_locus
  path ('resmarker_microhaplotype_table.txt'), emit: microhaps
  path ('all_mutations_table.txt'), emit: new_mutations

  script:
  """
  python3 ${projectDir}/bin/resistance_marker_module.py \
    --allele_data_path ${alleledata} \
    --aligned_asv_table_path ${alignment_data} \
    --res_markers_info_path ${resmarkers} \
    --refseq_path ${refseq} \
    --n-cores ${task.cpus}
  """
}