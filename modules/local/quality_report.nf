process QUALITY_REPORT {
  
  label 'process_low'
  conda 'envs/qc-env.yml'

  publishDir(
      path: "${params.outDIR}",
      mode: 'copy',
      pattern: "{sample_coverage|amplicon_coverage}.txt"
  )

  publishDir(
      path: "${params.outDIR}/quality_report",
      mode: 'copy',
      pattern: 'quality_report/*'
  )

  input:
  path (sample_coverage) 
  path (amplicon_coverage)
  path (amplicon_info)
  
  output:
  file ('sample_coverage.txt')
  file ('amplicon_coverage.txt')
  file ('quality_report')

  shell:
  """

  # Rename input files to published versions
  test -f sample_coverage.txt || mv $sample_coverage sample_coverage.txt
  test -f amplicon_coverage.txt || mv $amplicon_coverage amplicon_coverage.txt

  test -d quality_report || mkdir quality_report
  Rscript ${projectDir}/bin/cutadapt_summaryplots.R amplicon_coverage.txt sample_coverage.txt $amplicon_info quality_report
  """
}