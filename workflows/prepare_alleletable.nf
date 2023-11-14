// modules
include { BUILD_ALLELETABLE } from '../modules/local/build_alleletable.nf'
include { CALCULATE_IDENTITY_SCORE } from '../modules/local/calculate_identity_score.nf'
include { BUILD_PSEUDOCIGAR } from '../modules/local/build_pseudocigar.nf'

workflow PREPARE_ALLELETABLE {

  take:
  denoise_ch
  alignments_ch

  main:

  // Build the pseudocigar string from the unmasked alignments
  BUILD_PSEUDOCIGAR(
    alignments_ch
  )

  // Calculate the Identity Score
  CALCULATE_IDENTITY_SCORE(
    alignments_ch,
    BUILD_PSEUDOCIGAR.out.pseudocigar
  )

  // Finally create the final allele table
  BUILD_ALLELETABLE(
    BUILD_PSEUDOCIGAR.out.pseudocigar,
    CALCULATE_IDENTITY_SCORE.out.identity_score
  )

  emit:
  alleledata = BUILD_ALLELETABLE.out.alleledata
}