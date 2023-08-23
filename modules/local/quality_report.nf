process QUALITY_REPORT {
  
  label 'process_low'

  publishDir(
      path: "${params.outDIR}",
      mode: 'copy'
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
  mv $sample_coverage sample_coverage.txt
  mv $amplicon_coverage amplicon_coverage.txt

  test -d quality_report || mkdir quality_report
  Rscript ${projectDir}/bin/cutadapt_summaryplots.R amplicon_coverage.txt sample_coverage.txt $amplicon_info quality_report
  """
}