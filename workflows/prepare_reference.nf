#!/usr/bin/env nextflow

include { CHECK_FOR_PANEL_DEFINITIONS } from '../modules/local/check_for_panel_definitions.nf'
include { PREPARE_REFERENCE_SEQUENCES } from '../subworkflows/local/prepare_reference_sequences.nf'

workflow PREPARE_REFERENCE {

  take:
  refseq_fasta

  main:
  // create the reference if the user has not provided one (but has a genome), otherwise use the user file
  def reference = refseq_fasta ?: PREPARE_REFERENCE_SEQUENCES().reference_ch

  // Make sure that we have a valid reference for the panel
  CHECK_FOR_PANEL_DEFINITIONS(params.amplicon_info, reference)

  emit:
  reference_ch = reference
}

