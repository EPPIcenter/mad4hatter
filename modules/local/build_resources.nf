process BUILD_AMPLICON_INFO {
    input:
    val pools
    val amplicon_info_paths
    val amplicon_info_output

    publishDir(
        path: "${params.outDIR}/panel_information",
        mode: 'copy'
    )

    output:
        path "${amplicon_info_output}", emit: amplicon_info

    script:
    """
    python3 ${projectDir}/bin/build_amplicon_info.py \
        --pools ${pools} \
        --amplicon_info_paths ${amplicon_info_paths} \
        --amplicon_info_output_path ${amplicon_info_output}
    """
}