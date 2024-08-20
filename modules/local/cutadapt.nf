/*
 * STEP - CUTADAPT
 * Prepare the primer files from the given amplicon_info file
 */


process CUTADAPT {

  tag "$pair_id"
  label 'process_low'
  conda 'envs/cutadapt-env.yml'

  input:
  file fwd_primers
  file rev_primers
  tuple val(pair_id), file(reads)
  val cutadapt_minlen
  val sequencer
  val allowed_errors

  output:
  path("*.SAMPLEsummary.txt"), emit: sample_summary
  path("*.AMPLICONsummary.txt"), emit: amplicon_summary
  path('demultiplexed_fastqs'), emit: demultiplexed_fastqs
  path('adapter_dimers/*'), emit: adapter_dimers
  path('no_adapter_dimers/*'), emit: no_adapter_dimers
  path('fastp_reports/*'), emit: fastp_reports

  publishDir(
    path: "${params.outDIR}/cutadapt",
    mode: 'copy',
    pattern: 'adapter_dimers/*'
  )

  publishDir(
    path: "${params.outDIR}/cutadapt",
    mode: 'copy',
    pattern: 'no_adapter_dimers/*'
  )

  publishDir(
    path: "${params.outDIR}/fastp_reports",
    mode: 'copy',
    pattern: 'fastp_reports/*'
  )

  script:
  """
  bash cutadapt_process.sh \
      -1 ${reads[0]} \
      -2 ${reads[1]} \
      -r ${rev_primers} \
      -f ${fwd_primers} \
      -m ${cutadapt_minlen} \
      -s ${sequencer} \
      -e ${allowed_errors} \
      -c ${task.cpus} \
      -o demultiplexed_fastqs
  """
}

