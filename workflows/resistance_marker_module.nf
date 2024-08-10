/*
 * WORKFLOW - RESISTANCE_MARKER_MODULE
 * 
 * This workflow is comprised of processes to build a table of resistance markers
 * of interest.
 */

include { BUILD_RESISTANCE_TABLE } from '../modules/local/build_resistance_table.nf'
include { BUILD_RESMARKER_INFO } from '../modules/local/build_resources.nf'

workflow RESISTANCE_MARKER_MODULE {

    take:
    // resmarkers_amplicon
    amplicon_info
    allele_data
    alignment_data
    reference

    main:
    
    def resmarkers_amplicon = null
    if ( params.resmarker_info == null ) {
        BUILD_RESMARKER_INFO(amplicon_info, params.principal_resmarkers, 'resmarker_info.tsv')
        resmarkers_amplicon = BUILD_RESMARKER_INFO.out.resmarker_info
    }
    else {
        resmarkers_amplicon = params.resmarkers_amplicon
    }

    BUILD_RESISTANCE_TABLE(
        allele_data,
        alignment_data,
        resmarkers_amplicon,
        reference
    )
}