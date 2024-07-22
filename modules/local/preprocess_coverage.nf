// modules/local/preprocess_coverage.nf

process PREPROCESS_COVERAGE {

  label 'process_single'

  input:
  path too_short_coverages
  path sample_coverages
  path amplicon_coverages
  
  output:
  path 'sample_coverage.txt', emit: sample_coverage
  path 'amplicon_coverage.txt', emit: amplicon_coverage
  path 'too_short.txt', emit: too_short_coverage

  script:
  """
  add_sample_name_column() {
    awk -v fname=\$(basename "\$1" | sed -e 's/.SAMPLEsummary.txt//g' -e 's/.AMPLICONsummary.txt//g') -v OFS="\\t" '{print fname, \$0}' "\$1"
  }
  
  echo -e "SampleID\\tStage\\tReads" > sample_coverage.txt
  echo -e "SampleID\\tLocus\\tReads" > amplicon_coverage.txt

  for file in \$(ls $sample_coverages)
  do
      add_sample_name_column \$file >> sample_coverage.txt
  done

  for file in \$(ls $amplicon_coverages)
  do
      add_sample_name_column \$file >> amplicon_coverage.txt
  done

  for file in \$(ls $too_short_coverages)
  do
      add_sample_name_column \$file >> too_short.txt
  done
  """
}
