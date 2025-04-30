process BUILD_AMPLICON_INFO {
    input:
    val pools
    val amplicon_info_paths
    val amplicon_info_output

    publishDir(
        path: "${params.outDIR}/panel_information",
        mode: 'copy'
    )
    publishDir(path: "${params.outDIR}/panel_information", mode: 'copy')
    publishDir(path: "${params.outDIR}/standard_outputs/", mode: 'copy')

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

process BUILD_TARGETED_REFERENCE {
    input:
    val reference_input_paths
    val reference_output_path

    publishDir(
        path: "${params.outDIR}/panel_information",
        mode: 'copy'
    )
    output:
        path "${reference_output_path}", emit: reference_fasta

    script:
    """
    python3  ${projectDir}/bin/merge_fasta.py \
        --reference_paths ${reference_input_paths} \
        --reference_output_path ${reference_output_path}
    """
}

process BUILD_RESMARKER_INFO {
    input:
        val amplicon_info
        path principal_resmarkers
        val resmarker_info_output_path

    publishDir(
        path: "${params.outDIR}/panel_information",
        mode: 'copy'
    )
    output:
        path "${resmarker_info_output_path}", emit: resmarker_info, optional: true

    script:
    """
    python3 ${projectDir}/bin/build_resmarker_info.py \
        --amplicon_info ${amplicon_info} \
        --principal_resmarkers ${principal_resmarkers} \
        --resmarker_info_output_path ${resmarker_info_output_path}
    """
}