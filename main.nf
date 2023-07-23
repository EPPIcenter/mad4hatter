#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

if ( params.readDIR == null ) {
  exit 0, "ERROR: readDIR must be specified."
}

if ( params.target == null ) {
  exit 0, "ERROR: target must be specified."
}

// Expand user directory if exists
outDIR = "${params.outDIR}".replaceFirst("^~", System.getProperty("user.home"))
readDIR = "${params.readDIR}".replaceFirst("^~", System.getProperty("user.home"))

// Set boilerplate parameters
params.QC_only         = false
params.reads           = "${readDIR}/*_R{1,2}*.fastq.gz"
params.amplicon_info   = "$projectDir/resources/${params.target}/${params.target}_amplicon_info.tsv"
params.scriptDIR       = "$projectDir/R_code"
// pf3D7_index            = "$projectDir/resources/${params.target}/3D7_ampseq"
// codontable             = "$projectDir/resources/${params.target}/codontable.txt"
params.resmarkers_amplicon    = "$projectDir/resources/${params.target}/resistance_markers_amplicon_${params.target}.txt"
params.codontable      = "$projectDir/templates/codontable.txt"

// Files

cutadapt_minlen = params.cutadapt_minlen

/*
Create 'read_pairs' channel that emits for each read pair a
tuple containing 3 elements: pair_id, R1, R2
*/

include { CREATE_PRIMER_FILES } from './modules/local/create_primer_files.nf'
include { CUTADAPT } from './modules/local/cutadapt.nf'
include { QUALITY_REPORT } from './modules/local/quality_report.nf'
include { DADA2_ANALYSIS } from './modules/local/dada2_analysis.nf'

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
    params.allowed_errors,
    params.cutadapt_cores
  )

  emit:
  sample_summary_ch = CUTADAPT.out.sample_summary
  amplicon_summary_ch = CUTADAPT.out.amplicon_summary
  demux_fastqs_ch = CUTADAPT.out.demultiplexed_fastqs
}

workflow QUALITY_CONTROL {
  
  take:
  sample_summary_ch
  amplicon_summary_ch

  main:
  QUALITY_REPORT(
    sample_summary_ch.collect(),
    amplicon_summary_ch.collect(),
    params.amplicon_info,
    params.qc_cores
  )
}

workflow DENOISE_AMPLICONS {

  take:
  demultiplexed_fastqs

  main:

  DADA2_ANALYSIS(
    demultiplexed_fastqs.collect(),
    params.amplicon_info,
    params.pool,
    params.band_size,
    params.omega_a,
    params.maxEE,
    params.just_concatenate
  )
}

workflow {
  read_pairs = channel.fromFilePairs( params.reads, checkIfExists: true )
  DEMULTIPLEX_AMPLICONS(read_pairs)
  if (params.QC_only == true) {

    // create a quality report with the raw data
    QUALITY_CONTROL(
      DEMULTIPLEX_AMPLICONS.out.sample_summary_ch,
      DEMULTIPLEX_AMPLICONS.out.amplicon_summary_ch
    )

    // exit here
    return
  }

  DENOISE_AMPLICONS(
    DEMULTIPLEX_AMPLICONS.out.demux_fastqs_ch
  )
}

