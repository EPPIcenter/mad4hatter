/*
 * WORKFLOW - RESISTANCE_MARKER_MODULE
 * 
 * This workflow is comprised of processes to build a table of resistance markers
 * of interest.
 */

include { BUILD_RESISTANCE_TABLE } from '../modules/local/build_resistance_table.nf'

workflow RESISTANCE_MARKER_MODULE {

    take:

    allele_data
    reference

    main:

    // Additional Modules (not part of the main workflow)
    BUILD_RESISTANCE_TABLE(
        allele_data,
        params.resmarkers_amplicon,
        reference
    )
}