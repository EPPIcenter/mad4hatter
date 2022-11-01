#!/usr/bin/env nextflow

if ( params.readDIR == null ) {
  exit 0, "ERROR: readDIR must be specified."
}

if ( params.target == null ) {
  exit 0, "ERROR: target must be specified."
}

// Expand user directory if exists
outDIR = "${params.outDIR}".replaceFirst("^~", System.getProperty("user.home"))
readDIR = "${params.readDIR}".replaceFirst("^~", System.getProperty("user.home"))

// Set boilerplate parameters
params.reads           = "${readDIR}/*_R{1,2}*.fastq.gz"
params.amplicon_info   = "$projectDir/resources/${params.target}/${params.target}_amplicon_info.tsv"
params.refseq_fasta    = "$projectDir/resources/${params.target}/${params.target}_refseq.fasta"
params.scriptDIR       = "$projectDir/R_code"

// Files
amplicon_info = file( params.amplicon_info )

QC_only = params.QC_only

cutadapt_minlen = params.cutadapt_minlen
if ( params.sequencer == 'miseq' ) { qualfilter = '--trim-n -q 10' } else { qualfilter = '--nextseq-trim=20' }

/*
Create 'read_pairs' channel that emits for each read pair a
tuple containing 3 elements: pair_id, R1, R2
*/
Channel
    .fromFilePairs( params.reads )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { read_pairs }

//Split the Channel
read_pairs.into { read_pairs_cutadapt; read_pairs_qualitycheck }

// cutadapt based quality filtering and QC
// Filter to only reads with primer detected and trim poly-g
process cutadapt {
        container "aarandad/ampseq_workflow:latest"

        conda (params.enable_conda ? 'envs/cutadapt-env.yml' : null)

        input:
        tuple pair_id, file(reads) from read_pairs_cutadapt
        file amplicon_info
        val cutadapt_minlen
        val qualfilter

        output:
        file('trimmed_demuxed') into dada2
        file '*{AMPLICON,SAMPLE}summary.txt' into cutadapt_report

        script:
        """
        #!/usr/bin/env bash
        set -e

        fwd_primers_file="fwd_primers.fasta"
        rev_primers_file="rev_primers.fasta"

        cat $amplicon_info | awk 'NR==1 {
          for (i = 1; i <= NF; i++) {
            if ( \$i == "amplicon" ) {amplicon=i}
            if ( \$i == "fwd_primer" ) {fwd_primer=i}
            if ( \$i == "rev_primer" ) {rev_primer=i}
          }
        } NR>=1 {
          print \$amplicon, \$fwd_primer, \$rev_primer
        }' > adapters.txt

        # sanity check. there should be 3 fields.
        nfields="\$(head -n 1 adapters.txt | awk '{print NF}')"
        if [[ \$nfields != 3 ]]; then echo "ERROR: Must have '\'fwd_primer\' and \'rev_primer\' in vX_amplicon_info file!!!"; exit 1; fi

        # note: assumes order (1) amplicon, (2) fwd_adpater, (3) rev_adapter

        echo 'NR>1 { printf(">%s\\n^%s\\n", \$1, \$2) }' > modify.awk
        cat adapters.txt | cut -d ' ' -f 1,2 > fwd.txt
        cat adapters.txt | cut -d ' ' -f 1,3 > rev.txt

        awk -f modify.awk fwd.txt > \$fwd_primers_file
        awk -f modify.awk rev.txt > \$rev_primers_file

        mkdir trimmed_demuxed
        mkdir trimmed_demuxed_unknown
        mkdir filtered_out
        mkdir filtered_in

        cutadapt \
            --action=trim \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -e 0 \
            --no-indels \
            --minimum-length ${cutadapt_minlen} \
            -o filtered_out/${pair_id}_dimer_R1.fastq.gz \
            -p filtered_out/${pair_id}_dimer_R2.fastq.gz \
            --untrimmed-output filtered_in/${pair_id}_filtered_R1.fastq.gz \
            --untrimmed-paired-output filtered_in/${pair_id}_filtered_R2.fastq.gz \
            ${reads[0]} \
            ${reads[1]} 1> ${pair_id}.cutadapt1_summary.txt

        grep -E "Total read pairs processed:" ${pair_id}.cutadapt1_summary.txt | paste -s -d";\n"  | sed 's/Total read pairs processed://' | sed 's/ //g' | sed 's/,//g' \
           | awk -v var1=${pair_id} 'BEGIN{FS=" ";OFS="\t"}{print var1, "Input",\$1};' > ${pair_id}.SAMPLEsummary.txt

        grep -E "Pairs discarded as untrimmed:" ${pair_id}.cutadapt1_summary.txt | paste -s -d";\n"  | sed 's/Pairs discarded as untrimmed://' | sed -r 's/[(].*//' | sed 's/ //g' | sed 's/,//g' \
           | awk -v var1=${pair_id} 'BEGIN{FS=" ";OFS="\t"}{print var1, "No Dimers",\$1};' >> ${pair_id}.SAMPLEsummary.txt

        cutadapt \
            --action=trim \
            -g file:\${fwd_primers_file} \
            -G file:\${rev_primers_file} \
            --pair-adapters \
            -e 0 \
            --no-indels \
            ${qualfilter} \
            --minimum-length ${cutadapt_minlen} \
            -o trimmed_demuxed/{name}_${pair_id}_trimmed_R1.fastq.gz \
            -p trimmed_demuxed/{name}_${pair_id}_trimmed_R2.fastq.gz \
            --untrimmed-output trimmed_demuxed_unknown/${pair_id}_unknown_R1.fastq.gz \
            --untrimmed-paired-output trimmed_demuxed_unknown/${pair_id}_unknown_R2.fastq.gz \
            filtered_in/${pair_id}_filtered_R1.fastq.gz \
            filtered_in/${pair_id}_filtered_R2.fastq.gz 1> ${pair_id}.cutadapt2_summary.txt

        grep -E "Pairs written" ${pair_id}.cutadapt2_summary.txt | cut -f 2- -d ":" | sed -r 's/[(].*//' | sed 's/ //g' | sed 's/,//g' | awk -v var1=${pair_id} 'BEGIN{FS=" ";OFS="\t"}{print var1, "Amplicons",\$1};' >> ${pair_id}.SAMPLEsummary.txt

        grep -E "Adapter|Sequence" ${pair_id}.cutadapt2_summary.txt | paste -s -d";\n" | sed 's/=== //' | cut -f 1,5 -d ";" | grep First | cut -f 4,7 -d " " \
           | awk -v var1=${pair_id} 'BEGIN{FS=" ";OFS="\t"}{print var1,\$1,\$2};' > ${pair_id}.trim.AMPLICONsummary.txt

        if [ "\$(ls -A trimmed_demuxed/)" ]; then
           for afile in trimmed_demuxed/*trimmed_R1.fastq.gz; do primer=`echo \$afile | cut -f 2- -d '/' | cut -f 1-3 -d '_'`;  samplename=${pair_id};  readcount=`zcat \$afile | awk 'NR % 4 ==1' | wc -l`; printf "%s\t%s\t%s\n" \$samplename \$primer \$readcount >> ${pair_id}.filt.AMPLICONsummary.txt; done
           else
           touch ${pair_id}.filt.AMPLICONsummary.txt
        fi
        """
}

// Quality checks on the cutadapt summary file
process qualitycheck {
        container "aarandad/ampseq_workflow:latest"

        conda (params.enable_conda ? 'pandoc' : null)

        publishDir "${outDIR}",
               saveAs: { filename -> "${filename}"
               }

        input:
        tuple pair_id, file(reads) from read_pairs_qualitycheck
        file summfile from cutadapt_report.collect().ifEmpty([])
        file amplicon_info

        output:
        file '*coverage.txt' into qualitycheck_report
        file('quality_report')

        
        script:
        """
        #!/usr/bin/env bash
        set -e

        echo $summfile | tr ' ' '\n' | grep 'filt.AMPLICONsummary.txt' | tr '\n' ' ' | xargs cat > ALL_filt.AMPLICONsummary.txt
        echo $summfile | tr ' ' '\n' | grep 'trim.AMPLICONsummary.txt' | tr '\n' ' ' | xargs cat > ALL_trim.AMPLICONsummary.txt
        echo $summfile | tr ' ' '\n' | grep 'SAMPLEsummary.txt' | tr '\n' ' ' | xargs cat > sample_coverage.txt

        if [ -s ALL_filt.AMPLICONsummary.txt ]; then
        awk 'NR == FNR { key[\$1,\$2] = \$3; next } { \$3 = ((\$1,\$2) in key) ? key[\$1,\$2] : 0 };1' OFS="\t"  ALL_filt.AMPLICONsummary.txt ALL_trim.AMPLICONsummary.txt > amplicon_coverage.txt
        else
        awk 'BEGIN{FS=OFS="\t"} { print \$1,\$2,0; }'  ALL_trim.AMPLICONsummary.txt > amplicon_coverage.txt
        fi

        mkdir quality_report
        Rscript ${params.scriptDIR}/cutadapt_summaryplots.R amplicon_coverage.txt sample_coverage.txt ${amplicon_info} quality_report
        """
}

// Dada2

process dada2_analysis {
        container "aarandad/ampseq_workflow:latest"

        publishDir "${outDIR}",
               saveAs: { filename -> "${filename}"
               }

        input:
        file 'trimmed_demuxed' from dada2.collect().ifEmpty([])
        file amplicon_info

        output:
        file '*.RData' into dada2_summary
        
        when : QC_only != "T"

        script:
        treat_no_overlap_differently = 'T'

        """
        Rscript ${params.scriptDIR}/dada_overlaps.R ${trimmed_demuxed} ${amplicon_info} $treat_no_overlap_differently dada2_output.RData
        """
}

// Dada2 Postprocessing
process dada2_postproc {

        container "aarandad/ampseq_workflow:latest"

        publishDir "${outDIR}",
               saveAs: { filename -> "${filename}"
               }

        input:
        file rdatafile from dada2_summary

        output:
        file '*.{RDS,txt,csv}' into dada2_proc
        
        when : QC_only != "T"

        script:

        """
        Rscript ${params.scriptDIR}/postdada_rearrange.R $rdatafile ${params.homopolymer_threshold} ${params.refseq_fasta}
        """
}

