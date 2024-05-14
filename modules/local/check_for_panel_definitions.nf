/*
 * STEP - CHECK_FOR_PANEL_DEFINITIONS
 * Ensure that are amplicons in the panel have a reference sequence 
 */

process CHECK_FOR_PANEL_DEFINITIONS {

  label 'process_single'
  errorStrategy 'terminate'

  input:
  path amplicon_info
  path refseq_fasta

  script:
  """
  Rscript ${projectDir}/bin/check_for_panel_definitions.R \
    --ampliconFILE ${amplicon_info} \
    --reference-fasta ${refseq_fasta} \
  """
}