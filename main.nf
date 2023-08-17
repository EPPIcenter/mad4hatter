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
params.resmarkers_amplicon    = "$projectDir/resources/${params.target}/resistance_markers_amplicon_${params.target}.txt"

// Files
cutadapt_minlen = params.cutadapt_minlen

/*
Create 'read_pairs' channel that emits for each read pair a
tuple containing 3 elements: pair_id, R1, R2
*/

// workflows
include { DEMULTIPLEX_AMPLICONS } from './workflows/demultiplex_amplicons.nf'
include { DENOISE_AMPLICONS_1 } from './workflows/denoise_amplicons_1.nf'
include { DENOISE_AMPLICONS_2 } from './workflows/denoise_amplicons_2.nf'
include { RESISTANCE_MARKER_MODULE } from './workflows/resistance_marker_module.nf'
include { QUALITY_CONTROL} from './workflows/quality_control.nf'


// modules
include { BUILD_ALLELETABLE } from './modules/local/build_alleletable.nf'

// main workflow
workflow {
  read_pairs = channel.fromFilePairs( params.reads, checkIfExists: true )
  DEMULTIPLEX_AMPLICONS(read_pairs)

  // if QC_only is set, generate the report and exit early
  if (params.QC_only == true) {

    // create a quality report with the raw data
    QUALITY_CONTROL(
      DEMULTIPLEX_AMPLICONS.out.sample_summary_ch,
      DEMULTIPLEX_AMPLICONS.out.amplicon_summary_ch
    )

    // exit here
    return
  }

  // Denoising (DADA specific)
  DENOISE_AMPLICONS_1(
    DEMULTIPLEX_AMPLICONS.out.demux_fastqs_ch
  )

  // Masking, collapsing ASVs
  DENOISE_AMPLICONS_2(
    DENOISE_AMPLICONS_1.out.denoise_ch
  )

  // Finally create the final allele table
  BUILD_ALLELETABLE(
    DENOISE_AMPLICONS_1.out.denoise_ch,
    DENOISE_AMPLICONS_2.out.results_ch,
  )

  // Create the quality report now
  QUALITY_CONTROL(
    DEMULTIPLEX_AMPLICONS.out.sample_summary_ch,
    DEMULTIPLEX_AMPLICONS.out.amplicon_summary_ch
  )

  // By default, run the resistance marker module in the main workflow
  RESISTANCE_MARKER_MODULE(
    BUILD_ALLELETABLE.out.alleledata,
    DENOISE_AMPLICONS_2.out.reference_ch
  )
}

