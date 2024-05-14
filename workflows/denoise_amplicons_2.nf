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

workflow DENOISE_AMPLICONS_2 {

  take: 
  denoise_ch
  reference_ch
  main:

  // custom trimming
  denoise_ch = params.just_concatenate ? 
    COLLAPSE_CONCATENATED_READS(denoise_ch) : denoise_ch

  // use the denoised sequences and align them to the reference
  ALIGN_TO_REFERENCE(
    denoise_ch,
    reference_ch,
    params.amplicon_info
  )

  FILTER_ASVS(
    ALIGN_TO_REFERENCE.out.alignments
  )

  if (params.masked_fasta == null && (params.mask_tandem_repeats || params.mask_homopolymers)) {
    MASK_LOW_COMPLEXITY_REGIONS(
      reference_ch,
      FILTER_ASVS.out.filtered_alignments_ch
    )

    BUILD_PSEUDOCIGAR(
      MASK_LOW_COMPLEXITY_REGIONS.out.masked_alignments
    )
  } else {
    // Build the pseudocigar string from the unmasked alignments
    BUILD_PSEUDOCIGAR(
      FILTER_ASVS.out.filtered_alignments_ch
    )
  }

  emit:
  results_ch = BUILD_PSEUDOCIGAR.out.pseudocigar
}