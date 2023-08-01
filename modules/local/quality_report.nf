process QUALITY_REPORT {
  
  publishDir(
      path: "${params.outDIR}",
      mode: 'copy'
  )

  input:
  path (sample_coverage) 
  path (amplicon_coverage)
  path (amplicon_info)
  val (cores)
  
  output:
  file ('sample_coverage.txt')
  file ('amplicon_coverage.txt')
  file ('quality_report')

  shell:

  """
  # create a function to add the sample name as a column
  add_sample_name_column() {
    awk -v fname=\$(basename "\$1" | sed -e 's/.SAMPLEsummary.txt//g' -e 's/.AMPLICONsummary.txt//g') -v OFS="\\t" '{print fname, \$0}' "\$1"
  }

  export -f add_sample_name_column

  # add headers for easier readability
  echo -e "SampleName\\tX\\tNumReads" > sample_coverage.txt
  echo -e "SampleName\\tAmplicon\\tNumReads" > amplicon_coverage.txt

  # concatenate sample and amplicon coverage files form all samples
  printf "%s\n" ${sample_coverage.join(' ')} | parallel --jobs $cores add_sample_name_column {} >> sample_coverage.txt
  printf "%s\n" ${amplicon_coverage.join(' ')} | parallel --jobs $cores add_sample_name_column {} >> amplicon_coverage.txt

  test -d quality_report || mkdir quality_report
  Rscript ${projectDir}/bin/cutadapt_summaryplots.R amplicon_coverage.txt sample_coverage.txt ${amplicon_info} quality_report
  """
}