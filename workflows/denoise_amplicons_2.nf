
include { ALIGN_TO_REFERENCE } from '../modules/local/align_to_reference.nf'
include { MASK_LOW_COMPLEXITY_REGIONS } from '../subworkflows/local/mask_low_complexity_regions.nf'
include { PREPARE_REFERENCE_SEQUENCES } from '../subworkflows/local/prepare_reference_sequences.nf'


workflow DENOISE_AMPLICONS_2 {

  take: 
  denoise_ch

  main:

  // create the reference if the user has not provided one (but has a genome), otherwise use the user file
  def reference = (params.refseq_fasta == null) ? PREPARE_REFERENCE_SEQUENCES().reference_ch : params.refseq_fasta

  println(reference)

  // use the denoised sequences and align them to the reference
  ALIGN_TO_REFERENCE(
    denoise_ch,
    reference,
    params.amplicon_info,
    params.parallel,
    params.n_cores
  )

  if (params.masked_fasta == null && (params.mask_tandem_repeats || params.mask_homopolymers)) {
    MASK_LOW_COMPLEXITY_REGIONS(
      reference,
      ALIGN_TO_REFERENCE.out.alignments
    )
  }

  // BUILD_PSEUDO_CIGAR

  // emit:
  // masked_allele_data_ch = MASK_SEQUENCES.out.allele_data
}