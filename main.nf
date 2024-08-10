#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Expand user directory if exists
outDIR = "${params.outDIR}".replaceFirst("^~", System.getProperty("user.home"))
readDIR = "${params.readDIR}".replaceFirst("^~", System.getProperty("user.home"))

// Set boilerplate parameters
params.QC_only         = false
params.reads           = "${readDIR}/*_R{1,2}*.fastq.gz"
params.help	       = false

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
include { GENERATE_AMPLICON_INFO } from './workflows/process_inputs.nf'
include { VALIDATE_INPUTS } from './workflows/validate_inputs.nf'
include { BUILD_RESMARKER_INFO } from './modules/local/build_resources.nf'

// workflows
include { QC_ONLY } from './workflows/qc_only.nf'
include { POSTPROC_ONLY } from './workflows/postproc_only.nf'


// modules
include { BUILD_ALLELETABLE } from './modules/local/build_alleletable.nf'

def helpMessage() {
  log.info """
    Usage:
      The typical command for running the pipeline is as follows:
      nextflow run main.nf --readDIR data/testdata --target v4

    Mandatory arguments:
      --pools     The pools that were used for sequencing. [Options: 1A,1B,2,5]
      --sequencer	The sequencer used to produce your data. [Options: miseq, nextseq]

    Mandatory for `complete` (default) or `qc` workflow: 
      --readDIR		Path to folder containing fastq files

    Mandatory for `postprocessing` workflow: 
      --denoised_asvs           Path to denoised ASVs from DADA2. Used to only run the postprocessing workflow

    Optional arguments:
      --outDIR                  Path to folder to place final output [Default: results]
      --workflow                Workflow option to be run [Options: complete (default), qc, postprocessing]

      (Nextflow parameters. Note the flags have "-" and not "--" at the start) 
      -profile                  Runtime profile [Options: sge,apptainer, docker, conda]
      -config                   Resource configurations for each process

      (DADA2 parameters)
      --omega_a                 Level of statistical evidence required for DADA2 to infer a new ASV [Default: 1e-120]
      --pool                    Pooling method for DADA2 to process ASVs [Options: pseudo (default), true, false]
      --band_size               Limit on the net cumulative number of insertions of one sequence relative to the other in DADA2 [Default: 16]
      --maxEE                   Limit on number of expected errors within a read during filtering and trimming within DADA2 [Default: 3]

      (Post processing parameters)
      --concat_non_overlaps     Whether to concatenate or discard any sequences that DADA2 was unable to be merge 
      --refseq_fasta            Path to targeted reference 
      --genome                  Path to full genome covering all targets 
      --homopolymer_threshold   The length a homopolymer must reach to be masked [Default: 5]
      --trf_min_score           The alignment of a tandem repeat must meet or exceed this alignment score to be masked [Default: 25]
      --trf_max_period          The pattern size must be less than this to be masked [Default: 3]

    Examples:
      More advanced usage 
      nnextflow run main.nf --readDIR data/testdata --outDIR results --pools 1A,1B,2 --sequencer nextseq -config custom.config

      Change runtime param to use Docker
      nextflow run main.nf --readDIR data/testdata --outDIR results --pools 1A,1B,2 -profile docker
      
      Only run the QC workflow
      nextflow run main.nf --readDIR data/testdata --outDIR results --pools 1A,1B,2 --workflow qc
      
      Only run the Postprocessing workflow
      nextflow run main.nf --workflow postprocessing --pools 1A,1B,2 --denoised_asvs results/raw_dada2_output/dada2.clusters.txt --outDIR postprocessing_results
      
      Alter Dada2 params 
      nextflow run main.nf --readDIR data/testdata --outDIR results --pools 1A,1B,2 --omega_a 1e-40 --pool false --band_size 20 --maxEE 4
      
      Set full genome 
      nextflow run main.nf --readDIR data/testdata --outDIR results --pools 1A,1B,2 --genome data/reference/v1/PkPfPmPoPv.fasta
      
      Set targeted reference 
      nextflow run main.nf --readDIR data/testdata --outDIR results --pools 1A,1B,2 --refseq_fasta resources/v4/v4_reference.fasta
      
      Alter Masking parameters 
      nextflow run main.nf --readDIR data/testdata --outDIR results --pools 1A,1B,2 --homopolymer_threshold 2 --trf_min_score 30 --trf_max_period 5
        """.stripIndent()
}

// main workflow
workflow {
  // Print help if requested
  if (params.help) {
      helpMessage()
      exit 0
  }

  VALIDATE_INPUTS()

  def amplicon_info = (params.amplicon_info == null) ? GENERATE_AMPLICON_INFO().amplicon_info_ch : params.amplicon_info
  // Convert workflow to lowercase
  def workflow = params.workflow?.toLowerCase()

  if (workflow=='qc') {
    // Run QC Only Workflow
    QC_ONLY(amplicon_info, params.reads)

  } else if (workflow=='postprocessing') {
    // Run Postprocessing only
    POSTPROC_ONLY(amplicon_info, params.denoised_asvs)

  } else if (workflow=='complete'){
    // Create read pairs channel from fastq data
    read_pairs = channel.fromFilePairs( params.reads, checkIfExists: true )

    // Trim and demultiplex amplicons by amplicon
    DEMULTIPLEX_AMPLICONS(amplicon_info, read_pairs)

    // Denoising (DADA specific)
    DENOISE_AMPLICONS_1(
      amplicon_info,
      DEMULTIPLEX_AMPLICONS.out.demux_fastqs_ch
    )

    // Masking, collapsing ASVs
    DENOISE_AMPLICONS_2(
      amplicon_info,
      DENOISE_AMPLICONS_1.out.denoise_ch
    )

    // Finally create the final allele table
    BUILD_ALLELETABLE(
      DENOISE_AMPLICONS_1.out.denoise_ch,
      DENOISE_AMPLICONS_2.out.results_ch
    )

    // Create the quality report now
    QUALITY_CONTROL(
      amplicon_info,
      DEMULTIPLEX_AMPLICONS.out.sample_summary_ch,
      DEMULTIPLEX_AMPLICONS.out.amplicon_summary_ch,
      BUILD_ALLELETABLE.out.alleledata,
      DENOISE_AMPLICONS_1.out.denoise_ch
    )
    // // RESMARKER
    // def resmarkers_amplicon = null
    // if ( params.resmarker_info == null ) {
    //     BUILD_RESMARKER_INFO(amplicon_info, params.principal_resmarkers, 'resmarker_info.tsv')
    //     resmarkers_amplicon = BUILD_RESMARKER_INFO.out.resmarker_info
    // }
    // else {
    //     resmarkers_amplicon = params.resmarkers_amplicon
    // }

    RESISTANCE_MARKER_MODULE(
      // resmarkers_amplicon,
      amplicon_info,
      BUILD_ALLELETABLE.out.alleledata,
      DENOISE_AMPLICONS_2.out.aligned_asv_table,
      DENOISE_AMPLICONS_2.out.reference_ch
    )
  } else {
    exit 0, log.error("`--workflow` invalid.")
  }
}

workflow.onComplete {
  def outputDir = new File("${params.outDIR}/run")
  if (!outputDir.exists()) {
      outputDir.mkdirs()
  }

  record_params()
  record_runtime()
}


/* Record parmameters
 *
 * Records set parameters for the run
 *
 */
def record_params() {
    
    def output = new File("${params.outDIR}/run/parameters.tsv")

    params.each{ k, v -> 
        output.append("${k}\t${v}\n")
    }
}


/* Record runtime information
 *
 * Records runtime and environment information and writes summary to a tabulated (tsv) file
 *
 */
def record_runtime() {

    def output = new File("${params.outDIR}/run/runtime.tsv")

    // Append ContainerEngine first as shown in your example
    output.append("PipelineVersion\t${workflow.manifest.version}\n")
    output.append("ContainerEngine\t${workflow.containerEngine}\n")
    output.append("Duration\t${workflow.duration}\n")
    output.append("CommandLine\t${workflow.commandLine}\n")
    output.append("CommitId\t${workflow.commitId}\n")
    output.append("Complete\t${workflow.complete}\n")
    output.append("ConfigFiles\t${workflow.configFiles.join(', ')}\n")
    output.append("Container\t${workflow.container}\n")
    output.append("ErrorMessage\t${workflow.errorMessage}\n")
    output.append("ErrorReport\t${workflow.errorReport}\n")
    output.append("ExitStatus\t${workflow.exitStatus}\n")
    output.append("HomeDir\t${workflow.homeDir}\n")
    output.append("LaunchDir\t${workflow.launchDir}\n")
    output.append("Manifest\t${workflow.manifest}\n")
    output.append("Profile\t${workflow.profile}\n")
    output.append("ProjectDir\t${workflow.projectDir}\n")
    output.append("Repository\t${workflow.repository}\n")
    output.append("Resume\t${workflow.resume}\n")
    output.append("Revision\t${workflow.revision}\n")
    output.append("RunName\t${workflow.runName}\n")
    output.append("ScriptFile\t${workflow.scriptFile}\n")
    output.append("ScriptId\t${workflow.scriptId}\n")
    output.append("ScriptName\t${workflow.scriptName}\n")
    output.append("SessionId\t${workflow.sessionId}\n")
    output.append("Start\t${workflow.start}\n")
    output.append("StubRun\t${workflow.stubRun}\n")
    output.append("Success\t${workflow.success}\n")
    output.append("UserName\t${workflow.userName}\n")
    output.append("WorkDir\t${workflow.workDir}\n")
    output.append("NextflowBuild\t${nextflow.build}\n")
    output.append("NextflowTimestamp\t${nextflow.timestamp}\n")
    output.append("NextflowVersion\t${nextflow.version}\n")
}
