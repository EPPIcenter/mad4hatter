
include { QUALITY_REPORT } from '../modules/local/quality_report.nf'

workflow QUALITY_CONTROL {
  
  take:
  sample_summary_ch
  amplicon_summary_ch

  main:
  QUALITY_REPORT(
    sample_summary_ch.collect(),
    amplicon_summary_ch.collect(),
    params.amplicon_info
  )
}