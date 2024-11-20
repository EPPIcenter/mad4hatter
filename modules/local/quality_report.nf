process QUALITY_REPORT {
  
  label 'process_low'
  conda 'envs/qc-env.yml'

  publishDir(
      path: "${params.outDIR}",
      mode: 'copy',
      pattern: '*_coverage.txt' // Matches both sample_coverage.txt and amplicon_coverage.txt
  )

  publishDir(
      path: "${params.outDIR}/quality_report",
      mode: 'copy'
  )

  input:
  path (sample_coverage) 
  path (amplicon_coverage)
  path (amplicon_info)
  
  output:
  path ('sample_coverage.txt')
  path ('amplicon_coverage.txt')
  path ('amplicon_stats.txt')
  path ('length_vs_reads.pdf')
  path ('QCplots.html')
  path ('QCplots.Rmd')
  path ('reads_histograms.pdf')
  path ('swarm_plots.pdf')

  shell:
  """

  # Rename input files to published versions
  test -f sample_coverage.txt || mv $sample_coverage sample_coverage.txt
  test -f amplicon_coverage.txt || mv $amplicon_coverage amplicon_coverage.txt

  // test -d quality_report || mkdir quality_report
  Rscript ${projectDir}/bin/cutadapt_summaryplots.R amplicon_coverage.txt sample_coverage.txt $amplicon_info .
  """
}