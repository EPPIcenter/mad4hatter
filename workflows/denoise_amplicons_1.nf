/*
 * WORKFLOW - DENOISE_AMPLICONS_1
 * 
 * This workflow uses DADA2 to denoise amplicons
 */


include { DADA2_ANALYSIS } from '../modules/local/dada2_analysis.nf'

workflow DENOISE_AMPLICONS_1 {

  take:
  amplicon_info
  demultiplexed_fastqs

  main:

  // This module takes all amplicon-demultiplexed fastqs and runs DADA2
  DADA2_ANALYSIS(
    demultiplexed_fastqs.collect(),
    amplicon_info,
    params.dada2_pool,
    params.band_size,
    params.omega_a,
    params.maxEE,
    params.just_concatenate, 
    params.matchIDs
  )

  emit: 
  denoise_ch = DADA2_ANALYSIS.out.dada2_clusters
}

