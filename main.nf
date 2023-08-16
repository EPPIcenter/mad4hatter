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
include { IDENTIFY_VARIANTS } from './workflows/identify_variants.nf'
include { RESISTANCE_MARKER_MODULE } from './workflows/resistance_marker_module.nf'
include { QUALITY_CONTROL} from './workflows/quality_control.nf'
include { ALIGN_TO_REFERENCE } from './modules/local/align_to_reference.nf'
include { ASV_ALIGNMENT } from './workflows/asv_alignment.nf'

// modules
include { BUILD_ALLELETABLE } from './modules/local/build_alleletable.nf'

// main workflow
workflow {
  read_pairs = channel.fromFilePairs( params.reads, checkIfExists: true )
  DEMULTIPLEX_AMPLICONS(read_pairs)
  DEMULTIPLEX_AMPLICONS.out.demux_coverage_ch.view()

}


  // // if QC_only is set, generate the report and exit early
  // if (params.QC_only == true) {

  //   // create a quality report with the raw data
  //   QUALITY_CONTROL(
  //     DEMULTIPLEX_AMPLICONS.out.demux_coverage_ch
  //   )

  //   // exit here
  //   return
  // }


  // Denoising (DADA specific)
  // DENOISE_AMPLICONS_1(
  //   DEMULTIPLEX_AMPLICONS.out.demux_fastqs_ch
  // )

  

  // ASV_ALIGNMENT(
  //   DENOISE_AMPLICONS_1.out.denoise_ch
  // )

  // // Masking, collapsing ASVs
  // IDENTIFY_VARIANTS(
  //   ASV_ALIGNMENT.out.results_ch
  //   ASV_ALIGNMENT.OUT.
  // )

  // // Finally create the final allele table
  // BUILD_ALLELETABLE(
  //   DENOISE_AMPLICONS_1.out.denoise_ch,
  //   IDENTIFY_VARIANTS.out.results_ch,
  // )

  // // Create a merged channel for each coverage type
  // sample_coverage_combined = Channel
  //     .from([
  //         ['X', 'demux', DEMULTIPLEX_AMPLICONS.out.sample_summary_ch],
  //         ['X', 'denoised', COMPUTE_ASV_COVERAGE.out.sample_asv_coverage]
  //     ])
  //     // Flattening the data so the sample_coverage_ch will emit tuples like: ['X', 'demux', file_path]
  //     .flatMap { it -> it }

  // amplicon_coverage_combined = Channel
  //     .from([
  //         ['demux', DEMULTIPLEX_AMPLICONS.out.amplicon_summary_ch],
  //         ['denoised', COMPUTE_ASV_COVERAGE.out.amplicon_asv_coverage]
  //     ])

  // // Combine the sample coverage channels into a merged channel
  // sample_coverage_ch = sample_coverage_combined
  //     .mix()
  //     .groupTuple()
  //     .map { type, files -> [type, files.flatten()] }

  // // Combine the amplicon coverage channels into a merged channel
  // amplicon_coverage_ch = amplicon_coverage_combined
  //     .mix()
  //     .groupTuple()
  //     .map { type, files -> [type, files.flatten()] }


  // sample_coverage_ch.View()



//   // Create the quality report now
//   QUALITY_CONTROL(
//     DEMULTIPLEX_AMPLICONS.out.sample_summary_ch,
//     DEMULTIPLEX_AMPLICONS.out.amplicon_summary_ch,
//     ASV_ALIGNMENT.out.sample_asv_coverage,
//     ASV_ALIGNMENT.out.amplicon_asv_coverage
//   )

//   // By default, run the resistance marker module in the main workflow
//   RESISTANCE_MARKER_MODULE(
//     BUILD_ALLELETABLE.out.alleledata,
//     IDENTIFY_VARIANTS.out.reference_ch
//   )
// }

