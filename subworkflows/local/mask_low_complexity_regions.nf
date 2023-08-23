/*
 * WORKFLOW - ALIGN-TO-REFERENCE
 * 
 * This workflow allows users to optionally mask homopolymers and / or tandem repeats
 */


include { MASK_REFERENCE_TANDEM_REPEATS } from '../../modules/local/mask_reference_tandem_repeats.nf'
include { MASK_REFERENCE_HOMOPOLYMERS } from '../../modules/local/mask_reference_homopolymers.nf'
include { MASK_SEQUENCES } from '../../modules/local/mask_sequences.nf'

workflow MASK_LOW_COMPLEXITY_REGIONS {
  
  take:

  reference
  alignments  

  main:
  
  // optionally mask sequences
  if (params.mask_tandem_repeats) {
    MASK_REFERENCE_TANDEM_REPEATS(
      reference,
      params.trf_min_score,
      params.trf_max_period
    )

    mask_reference_tandem_repeats_ch = MASK_REFERENCE_TANDEM_REPEATS.out.masked_fasta
  }


  if (params.mask_homopolymers) {
    MASK_REFERENCE_HOMOPOLYMERS(
      reference,
      params.homopolymer_threshold
    )
    
    mask_reference_homopolymers_ch = MASK_REFERENCE_HOMOPOLYMERS.out.masked_fasta
  }

  masked_reference_ch = mask_reference_tandem_repeats_ch.mix(mask_reference_homopolymers_ch)

  MASK_SEQUENCES(
    masked_reference_ch.collect(),
    alignments
  )

  emit:

  masked_alignments = MASK_SEQUENCES.out.masked_alignments
}