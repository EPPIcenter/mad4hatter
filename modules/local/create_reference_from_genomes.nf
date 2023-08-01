/*
 * STEP - CREATE_REFERENCE_FROM_GENOMES
 * 
 */


process CREATE_REFERENCE_FROM_GENOMES {

    publishDir(
      path: "${params.outDIR}",
      mode: 'copy'
  )

  input:
  path genome
  path amplicon_info
  val refseq_fasta

  output:
  path "${refseq_fasta}", emit: reference_fasta

  script:
  """
  Rscript ${projectDir}/bin/create_reference_from_genomes.R --ampliconFILE ${amplicon_info} --genome ${genome} --output ${refseq_fasta}
  """
}