include { CREATE_PRIMER_FILES } from '../modules/local/create_primer_files.nf'
include { CUTADAPT } from '../modules/local/cutadapt.nf'
include { QUALITY_REPORT } from '../modules/local/quality_report.nf'

workflow DEMULTIPLEX_AMPLICONS {

  take: 
  read_pairs

  main:
  CREATE_PRIMER_FILES(params.amplicon_info)
  CUTADAPT(
    CREATE_PRIMER_FILES.out.fwd_primers,
    CREATE_PRIMER_FILES.out.rev_primers,
    read_pairs,
    params.cutadapt_minlen,
    params.sequencer,
    params.allowed_errors
  )

  emit:
  demux_fastqs_ch = CUTADAPT.out.demultiplexed_fastqs

  // Combine the coverage channels into a single channel
  demux_coverage_ch = Channel.from([
    ['demux_sample', CUTADAPT.out.sample_summary],
    ['demux_amplicon', CUTADAPT.out.amplicon_summary]
  ])
}

