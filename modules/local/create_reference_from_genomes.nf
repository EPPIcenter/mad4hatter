/*
 * STEP - CREATE_REFERENCE_FROM_GENOMES
 * 
 */


process CREATE_REFERENCE_FROM_GENOMES {

  label 'process_medium'
  conda 'envs/postproc-env.yml'

  input:
  path genome
  path amplicon_info
  val refseq_fasta

  publishDir(
        path: "${params.outDIR}/panel_information",
        mode: 'copy'
  )

  output:
  path "${refseq_fasta}", emit: reference_fasta

  script:
  """
  Rscript ${projectDir}/bin/create_reference_from_genomes.R \
    --ampliconFILE ${amplicon_info} \
    --genome ${genome} \
    --output ${refseq_fasta} \
    --ncores ${task.cpus}
  """
}