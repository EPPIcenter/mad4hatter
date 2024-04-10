// modules/local/postprocess_coverage.nf

process POSTPROCESS_COVERAGE {

  label 'process_single'

  input:
  path alleledata
  path clusters
  path sample_coverage
  path amplicon_coverage
  
  output:
  path 'sample_coverage_postprocessed.txt', emit: postprocess_sample_coverage
  path 'amplicon_coverage_postprocessed.txt', emit: postprocess_amplicon_coverage

  script:
  """
  # Ensure only one file is provided for each coverage type
  if [ \$(ls -l $sample_coverage | wc -l) -ne 1 ] || [ \$(ls -l $amplicon_coverage | wc -l) -ne 1 ]; then
    echo "Error: Multiple coverage files detected. Expecting only one file for each type."
    exit 1
  fi
  
  Rscript ${projectDir}/bin/asv_coverage.R \
    --alleledata $alleledata \
    --clusters $clusters \
    --sample-coverage $sample_coverage \
    --amplicon-coverage $amplicon_coverage
  """
}
