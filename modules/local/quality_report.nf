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
  # create a function to add the sample name as a column
  add_sample_name_column() {
    awk -v fname=\$(basename "\$1" | sed -e 's/.SAMPLEsummary.txt//g' -e 's/.AMPLICONsummary.txt//g') -v OFS="\\t" '{print fname, \$0}' "\$1"
  }

  export -f add_sample_name_column

  # add headers for easier readability
  echo -e "SampleID\\tX\\tReads" > sample_coverage.txt
  echo -e "SampleID\\tLocus\\tReads" > amplicon_coverage.txt

  # concatenate sample and amplicon coverage files form all samples
  for file in ${sample_coverage.join(' ')}
  do
    add_sample_name_column \$file >> sample_coverage.txt
  done

  for file in ${amplicon_coverage.join(' ')}
  do
    add_sample_name_column \$file >> amplicon_coverage.txt
  done

  test -d quality_report || mkdir quality_report
  Rscript ${projectDir}/bin/cutadapt_summaryplots.R amplicon_coverage.txt sample_coverage.txt ${amplicon_info} quality_report
  """
}