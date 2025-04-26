/*
 * WORKFLOW - DENOISE_AMPLICONS_2
 * 
 * This workflow uses is comprised of multiple postprocessing steps to reduce noise with masking, 
 * and identify difference in the ASVs given a reference
 */


include { ALIGN_TO_REFERENCE } from '../modules/local/align_to_reference.nf'
include { MASK_LOW_COMPLEXITY_REGIONS } from '../subworkflows/local/mask_low_complexity_regions.nf'
include { PREPARE_REFERENCE_SEQUENCES } from '../subworkflows/local/prepare_reference_sequences.nf'
include { BUILD_PSEUDOCIGAR } from '../modules/local/build_pseudocigar.nf'
include { FILTER_ASVS } from '../modules/local/filter_asvs.nf'
include { COLLAPSE_CONCATENATED_READS } from '../modules/local/collapse_concatenated_reads.nf'

workflow MASK_AMPLICONS {

  take: 
  reference
  alignments

  main:

  MASK_LOW_COMPLEXITY_REGIONS(
    reference,
    alignments
  )
  masked_alignment_table_ch = MASK_LOW_COMPLEXITY_REGIONS.out.masked_alignments
  BUILD_PSEUDOCIGAR( 
    masked_alignment_table_ch
  )
  
  emit:
  masked_pseudocigar_ch = BUILD_PSEUDOCIGAR.out.pseudocigar
  masked_reference_ch = reference
  masked_aligned_asv_table = masked_alignment_table_ch
}