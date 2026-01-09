process PUBLISH_COVERAGE {
  
  label 'process_single'

  publishDir(
      path: "${params.outDIR}",
      mode: 'copy'
  )

  input:
  path (sample_coverage) 
  path (amplicon_coverage)
  
  output:
  path 'sample_coverage.txt', emit: sample_coverage
  path 'amplicon_coverage.txt', emit: amplicon_coverage

  script:
  """
  test -f sample_coverage.txt || mv $sample_coverage sample_coverage.txt
  test -f amplicon_coverage.txt || mv $amplicon_coverage amplicon_coverage.txt
  """
}
