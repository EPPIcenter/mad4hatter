
include { CREATE_REFERENCE_FROM_GENOMES } from '../../modules/local/create_reference_from_genomes.nf'

workflow PREPARE_REFERENCE_SEQUENCES {
  take: 
  amplicon_info
  
  main:
  if (params.genome == null) {
    exit 1, log.error("You must specify a genome file.")
  }

  def genome = null

  if (params.genomes.containsKey(params.genome)) {
    // Use a predefined genome
    genome = params.genomes[params.genome].fasta
  } else {
    // Read in genome file provided by the user 
    File genome_file = new File(params.genome).absoluteFile
    if(!genome_file.exists()) {
        exit 1, log.error("The specified denoised_asvs file '${params.denoised_asvs}' does not exist.")
    }

    // Add debugging steps as this is user input
    log.debug("Genome path: ${params.genome}")
    log.debug("Absolute path: ${genome_file.absolutePath}")
    log.debug("Does it exist? ${genome_file.exists()}")

    // Create the Nextflow Channel
    genome = Channel.fromPath(genome_file)   
  }

  CREATE_REFERENCE_FROM_GENOMES(
    genome,
    amplicon_info,
    "reference.fasta" // TODO: update this to include pool names 
  )

  emit:
  reference_ch = (params.refseq_fasta == null) ? 
    CREATE_REFERENCE_FROM_GENOMES.out.reference_fasta :
    params.refseq_fasta
}