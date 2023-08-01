/*
 * STEP - LONG-AMPLICONS
 * Denoise the demultiplexed amplicon fastqs 
 */

process CUTADAPT {

    input:
    file fwd_primers
    file rev_primers
    tuple val(pair_id), file(reads)
    val cutadapt_minlen
    val sequencer
    val allowed_errors
    val cores

    output:
    path("*.SAMPLEsummary.txt"), emit: sample_summary
    path("*.AMPLICONsummary.txt"), emit: amplicon_summary
    path('demuliplexed_fastqs'), emit: demultiplexed_fastqs

    script:
    """
    Rscript ${params.projectDir}/bin/postdada_rearrange.R \
      --homopolymer-threshold ${homopolymer_threshold} \
      --refseq-fasta ${refseq_fasta} \
      --masked-fasta ${masked_fasta} \
      --n-cores ${params.n_cores} \
      --parallel \
      --sample-coverage ${sample_coverage} \
      --amplicon-coverage ${amplicon_coverage} \
      --amplicon-table ${amplicon_info} \
      --clusters \${clusters_rds}
    """
}
