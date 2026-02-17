/*
 * WORKFLOW - DENOISE_AMPLICONS_2
 * 
 * This workflow uses is comprised of multiple postprocessing steps to reduce noise with masking, 
 * and identify difference in the ASVs given a reference
 */


include { ALIGN_TO_REFERENCE } from '../modules/local/align_to_reference.nf'
include { MASK_LOW_COMPLEXITY_REGIONS } from '../subworkflows/local/mask_low_complexity_regions.nf'
include { PREPARE_REFERENCE_SEQUENCES } from '../subworkflows/local/prepare_reference_sequences.nf'
include { BUILD_PSEUDOCIGAR; BUILD_PSEUDOCIGAR as BUILD_MASKED_PSEUDOCIGAR } from '../modules/local/build_pseudocigar.nf'
include { FILTER_ASVS } from '../modules/local/filter_asvs.nf'
include { COLLAPSE_CONCATENATED_READS } from '../modules/local/collapse_concatenated_reads.nf'

workflow DENOISE_AMPLICONS_2 {

  take: 
  amplicon_info
  denoise_ch

  main:

  // custom trimming
  denoise_ch = params.just_concatenate ? 
    COLLAPSE_CONCATENATED_READS(denoise_ch) : denoise_ch

  // create the reference if the user has not provided one (but has a genome), otherwise use the user file
  def reference = (params.refseq_fasta == null) ? PREPARE_REFERENCE_SEQUENCES(amplicon_info).reference_ch : params.refseq_fasta

  // use the denoised sequences and align them to the reference
  ALIGN_TO_REFERENCE(
    denoise_ch,
    reference,
    amplicon_info
  )

  FILTER_ASVS(
    ALIGN_TO_REFERENCE.out.alignments
  )
  alignment_table_ch = FILTER_ASVS.out.filtered_alignments_ch

  // Build the pseudocigar string for the unmasked alignments
  BUILD_PSEUDOCIGAR(
      alignment_table_ch, 
      ".unmasked"
    )

  // Apply masking to the alignments
  if (params.masked_fasta == null && (params.mask_tandem_repeats || params.mask_homopolymers)) {
    MASK_LOW_COMPLEXITY_REGIONS(
      reference,
      FILTER_ASVS.out.filtered_alignments_ch
    )
    alignment_table_ch = MASK_LOW_COMPLEXITY_REGIONS.out.masked_alignments
  } 

  // Build the pseudocigar string for the masked alignments
  BUILD_MASKED_PSEUDOCIGAR(
    alignment_table_ch,
    ".masked"
  )
  
  emit:
  denoise_ch = denoise_ch 
  unmasked_pseudocigar = BUILD_PSEUDOCIGAR.out.pseudocigar
  masked_pseudocigar = BUILD_MASKED_PSEUDOCIGAR.out.pseudocigar
  reference_ch = reference
  aligned_asv_table = alignment_table_ch
}