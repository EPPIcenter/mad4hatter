include { CREATE_PRIMER_FILES } from '../modules/local/create_primer_files.nf'
include { CUTADAPT } from '../modules/local/cutadapt.nf'

// Remove the sample number and lane added by Illumina sequencers (e.g. SAMPLE_S1_L001 -> SAMPLE)
// Must match SEQUENCER_SUFFIX in bin/build_pmo.py
def standardise_sample_name(sample_name) {
  return sample_name.replaceFirst(/_S\d+(_L\d{3})?$/, '')
}

workflow DEMULTIPLEX_AMPLICONS {

  take: 
  amplicon_info
  read_pairs

  main:
  if (params.standardise_sample_names) {
    read_pairs = read_pairs
      .map { pair_id, reads -> tuple(standardise_sample_name(pair_id), reads) }
      .toList()
      .flatMap { pairs ->
        def duplicated = pairs.countBy { it[0] }.findAll { it.value > 1 }.keySet()
        if (duplicated) {
          error "Multiple fastq pairs have the same sample name after removing sequencer suffixes: ${duplicated.join(', ')}. " +
            "Merge lanes before running the pipeline or set `--standardise_sample_names false`."
        }
        pairs
      }
  }

  CREATE_PRIMER_FILES(amplicon_info)
  CUTADAPT(
    CREATE_PRIMER_FILES.out.fwd_primers,
    CREATE_PRIMER_FILES.out.rev_primers,
    read_pairs,
    params.cutadapt_minlen,
    params.gtrim,
    params.quality_score,
    params.allowed_errors
  )

  emit:
  sample_summary_ch = CUTADAPT.out.sample_summary
  amplicon_summary_ch = CUTADAPT.out.amplicon_summary
  demux_fastqs_ch = CUTADAPT.out.demultiplexed_fastqs
}
