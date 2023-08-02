/*
 * STEP - BUILD_ALLELETABLE
 * 
 */

process BUILD_ALLELETABLE {

  label 'process_single'

  publishDir(
    path: "${params.outDIR}",
    mode: 'copy'
  )

  input:
  path denoised_asvs 
  path processed_asvs

  output:
  path("allele_data.txt"), emit: alleledata

  script:
  """
  Rscript ${projectDir}/bin/build_alleletable.R \
    --denoised-asvs ${denoised_asvs} \
    --processed-asvs ${processed_asvs}
  """
}