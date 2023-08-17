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


workflow DENOISE_AMPLICONS_2 {

  take: 
  denoise_ch

  main:

  // create the reference if the user has not provided one (but has a genome), otherwise use the user file
  def reference = (params.refseq_fasta == null) ? PREPARE_REFERENCE_SEQUENCES().reference_ch : params.refseq_fasta

  // use the denoised sequences and align them to the reference
  ALIGN_TO_REFERENCE(
    denoise_ch,
    reference,
    params.amplicon_info
  )

  if (params.masked_fasta == null && (params.mask_tandem_repeats || params.mask_homopolymers)) {
    MASK_LOW_COMPLEXITY_REGIONS(
      reference,
      ALIGN_TO_REFERENCE.out.alignments
    )

    BUILD_PSEUDOCIGAR(
      MASK_LOW_COMPLEXITY_REGIONS.out.masked_alignments
    )
  } else {
    // Build the pseudocigar string from the unmasked alignments
    BUILD_PSEUDOCIGAR(
      ALIGN_TO_REFERENCE.out.alignments
    )
  }

  emit:
  results_ch = BUILD_PSEUDOCIGAR.out.pseudocigar
  reference_ch = reference
}