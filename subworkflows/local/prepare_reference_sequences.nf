
include { CREATE_REFERENCE_FROM_GENOMES } from '../../modules/local/create_reference_from_genomes.nf'
include { CONCATENATE_TARGETED_REFERENCE } from './concatenate_targeted_reference.nf'

workflow PREPARE_REFERENCE_SEQUENCES {
  take: 
  amplicon_info
  
  main:
  if (params.genome) {
    def genome = null
    // Read in genome file provided by the user 
    File genome_file = new File(params.genome).absoluteFile
    if(!genome_file.exists()) {
        exit 1, log.error("The specified genome file '${params.genome}' does not exist.")
    }

    // Add debugging steps as this is user input
    log.debug("Genome path: ${params.genome}")
    log.debug("Absolute path: ${genome_file.absolutePath}")
    log.debug("Does it exist? ${genome_file.exists()}")

    // Create the Nextflow Channel
    genome = Channel.fromPath(genome_file)   

    CREATE_REFERENCE_FROM_GENOMES(
      genome,
      amplicon_info,
      "reference.fasta" // TODO: update this to include pool names 
    )
  } else {
    CONCATENATE_TARGETED_REFERENCE()
  }
  emit:
  reference_ch = (params.genome == null) ? 
    CONCATENATE_TARGETED_REFERENCE.out.reference_fasta :
    CREATE_REFERENCE_FROM_GENOMES.out.reference_fasta 
}