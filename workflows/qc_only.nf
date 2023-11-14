
/*
 * WORKFLOW - QC_ONLY
 * 
 * This workflow is meant to drive quality control of amplicon sequencing data
 */

include { DEMULTIPLEX_AMPLICONS } from './demultiplex_amplicons.nf'
include { QUALITY_CONTROL} from './quality_control.nf'


workflow QC_ONLY {
  take:
  reads

  main:

  read_pairs = channel.fromFilePairs( params.reads, checkIfExists: true )

  DEMULTIPLEX_AMPLICONS(read_pairs)

  // create a quality report with the raw data
  QUALITY_CONTROL(
    DEMULTIPLEX_AMPLICONS.out.sample_summary_ch,
    DEMULTIPLEX_AMPLICONS.out.amplicon_summary_ch,
    null,
    null
  )
}