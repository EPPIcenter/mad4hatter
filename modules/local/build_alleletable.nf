/*
 * STEP - BUILD_ALLELETABLE
 * 
 */

process BUILD_ALLELETABLE {

  label 'process_single'
  conda 'envs/postproc-env.yml'

  publishDir(path: "${params.outDIR}/standard_outputs", pattern: "microhaplotypes_info.tsv",mode: 'copy')
  publishDir(path: "${params.outDIR}", pattern: "allele_data.txt", mode: 'copy')

  input:
  path amplicon_info
  path denoised_asvs 
  path processed_asvs
  path processed_asvs_unmasked
  path aligned_asv_table

  output:
  path("allele_data.txt"), emit: alleledata
  path("microhaplotypes_info.tsv"), emit: mhap_data

  script:
  """
  Rscript ${projectDir}/bin/build_alleletable.R \
    --amplicon-info ${amplicon_info} \
    --denoised-asvs ${denoised_asvs} \
    --processed-asvs ${processed_asvs} \
    --processed-asvs-unmasked ${processed_asvs_unmasked} \
    --aligned-asvs ${aligned_asv_table}
  """
}