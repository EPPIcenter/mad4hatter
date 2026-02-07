include { CREATE_PRIMER_FILES } from '../modules/local/create_primer_files.nf'
include { CUTADAPT } from '../modules/local/cutadapt.nf'

workflow DEMULTIPLEX_AMPLICONS {

  take: 
  amplicon_info
  read_pairs

  main:
  CREATE_PRIMER_FILES(amplicon_info)
  CUTADAPT(
    CREATE_PRIMER_FILES.out.fwd_primers,
    CREATE_PRIMER_FILES.out.rev_primers,
    read_pairs,
    params.cutadapt_minlen,
    params.gtrim,
    params.allowed_errors
  )

  emit:
  sample_summary_ch = CUTADAPT.out.sample_summary
  amplicon_summary_ch = CUTADAPT.out.amplicon_summary
  demux_fastqs_ch = CUTADAPT.out.demultiplexed_fastqs
}
