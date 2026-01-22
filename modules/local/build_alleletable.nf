/*
 * STEP - BUILD_ALLELETABLE
 * 
 */

process BUILD_ALLELETABLE {

  label 'process_single'
  conda 'envs/postproc-env.yml'

  publishDir(path: "${params.outDIR}/standard_output", pattern: "microhaplotypes_info*.tsv",mode: 'copy')
  publishDir(path: "${params.outDIR}", pattern: "allele_data.txt", mode: 'copy')

  input:
  path amplicon_info
  path denoised_asvs 
  path masked_pseudocigar_table
  path unmasked_pseudocigar_table
  path masked_asv_table

  output:
  path("allele_data.txt"), emit: alleledata
  path("microhaplotypes_info.tsv"), emit: microhaplotypes_info
  path("microhaplotypes_info_collapsed.tsv"), emit: microhaplotypes_info_collapsed

  script:
  """
  Rscript ${projectDir}/bin/build_alleletable.R \
    --amplicon-info ${amplicon_info} \
    --denoised-asvs ${denoised_asvs} \
    --masked-pseudocigar-table ${masked_pseudocigar_table} \
    --unmasked-pseudocigar-table ${unmasked_pseudocigar_table} \
    --masked-asv-table ${masked_asv_table}
  """
}
