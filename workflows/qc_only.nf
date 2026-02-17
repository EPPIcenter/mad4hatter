
/*
 * WORKFLOW - QC_ONLY
 * 
 * This workflow is meant to drive quality control of amplicon sequencing data
 */

include { DEMULTIPLEX_AMPLICONS } from './demultiplex_amplicons.nf'
include { QUALITY_CONTROL} from './quality_control.nf'


workflow QC_ONLY {
  take:
  amplicon_info
  reads

  main:

  read_pairs = channel.fromFilePairs( reads, checkIfExists: true )

  DEMULTIPLEX_AMPLICONS(amplicon_info, read_pairs)

  // create a quality report with the raw data
  QUALITY_CONTROL(
    DEMULTIPLEX_AMPLICONS.out.sample_summary_ch,
    DEMULTIPLEX_AMPLICONS.out.amplicon_summary_ch,
    null,
    null
  )

  emit:
    sample_coverage_ch = QUALITY_CONTROL.out.sample_coverage
    amplicon_coverage_ch = QUALITY_CONTROL.out.amplicon_coverage
}