include { GET_SPIKEIN_COUNTS } from '../modules/local/get-spikein-counts.nf'
include { PLOT_SPIKEIN_COUNTS } from '../modules/local/plot-spikein-counts.nf'
include { SPIKEIN_TRIM } from '../modules/local/spikein-trim.nf'
include { CREATE_PRIMER_FILES } from '../modules/local/create_primer_files.nf'
include { ALIGN_TO_REFERENCE } from '../modules/local/align_to_reference.nf'
include { CREATE_SPIKEIN_TABLE } from '../modules/local/create_spikein_table.nf'

workflow SPIKEIN_ANALYSIS {

  take:
  unknown_fastqs_ch
  amplicon_coverage_ch

  main:

  CREATE_PRIMER_FILES(params.spikein_primers)

  expected_spikein_ch = channel.fromPath( params.expected_spikein, checkIfExists: true )
  spikein_info_ch = channel.fromPath( params.spikein_info, checkIfExists: true )

  // demutltiplex spikeins
  SPIKEIN_TRIM(
    CREATE_PRIMER_FILES.out.fwd_primers,
    CREATE_PRIMER_FILES.out.rev_primers,
    unknown_fastqs_ch
  )

  // get the spikein counts from the trimmed spikin fastqs
  GET_SPIKEIN_COUNTS(
    SPIKEIN_TRIM.out.spikeins_demultiplexed
  )

  // plot the spikein counts in a heatmap
  PLOT_SPIKEIN_COUNTS(
    GET_SPIKEIN_COUNTS.out.spikein_counts.collect(),
    expected_spikein_ch,
    spikein_info_ch,
    amplicon_coverage_ch,
  )

  // create a spikein table
  // CREATE_SPIKEIN_TABLE(
  //   SPIKEIN_TRIM.out.spikeins_demultiplexed.collect()
  // )
}