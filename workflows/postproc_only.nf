/*
 * WORKFLOW - POSTPROC_ONLY
 * 
 * This workflow uses is comprised of multiple postprocessing steps to reduce noise with masking, 
 * and identify difference in the ASVs given a reference.
 */

include { DENOISE_AMPLICONS_2 } from './denoise_amplicons_2.nf'
include { BUILD_ALLELETABLE } from '../modules/local/build_alleletable.nf'

workflow POSTPROC_ONLY {

    take:
    amplicon_info
    denoised_asvs

    main:
    // Read in denoised asv file provided by the user 
    File denoised_asvs_file = new File(denoised_asvs).absoluteFile
    if(!denoised_asvs_file.exists()) {
        exit 1, "The specified denoised_asvs file '${denoised_asvs}' does not exist."
    }

    // Add debugging steps as this is user input
    log.debug("Denoised ASVs path: ${denoised_asvs}")
    log.debug("Absolute path: ${denoised_asvs_file.absolutePath}")
    log.debug("Does it exist? ${denoised_asvs_file.exists()}")

    // Create the Nextflow Channel
    denoise_ch = Channel.fromPath(denoised_asvs_file)
    
    // Run postprocessing only workflow
    DENOISE_AMPLICONS_2(amplicon_info, denoise_ch)

    // Allele table creation
    BUILD_ALLELETABLE(
      amplicon_info,
      denoise_ch,
      DENOISE_AMPLICONS_2.out.results_ch
    )
}