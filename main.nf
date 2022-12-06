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
params.QC_only         = false
params.reads           = "${readDIR}/*_R{1,2}*.fastq.gz"
params.amplicon_info   = "$projectDir/resources/${params.target}/${params.target}_amplicon_info.tsv"
params.scriptDIR       = "$projectDir/R_code"
// pf3D7_index            = "$projectDir/resources/${params.target}/3D7_ampseq"
// codontable             = "$projectDir/resources/${params.target}/codontable.txt"
params.resmarkers_amplicon    = "$projectDir/resources/${params.target}/resistance_markers_amplicon_${params.target}.txt"

// Files

cutadapt_minlen = params.cutadapt_minlen
if ( params.sequencer == 'miseq' ) { qualfilter = '--trim-n -q 10' } else { qualfilter = '--nextseq-trim=20' }

/*
Create 'read_pairs' channel that emits for each read pair a
tuple containing 3 elements: pair_id, R1, R2
*/

workflow {
  // declare variables upfront
  read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )

  // if this is null, a template codontable will be used
  codontable = params.codontable != null ? params.codontable : "$projectDir/templates/codontable.txt"


  // All workflows will produce a QC report (cutadapt needed by QC for length statistics)
  CUTADAPT(read_pairs_ch, params.amplicon_info, params.cutadapt_minlen, qualfilter)
  QUALITY_CHECK(read_pairs_ch, CUTADAPT.out[1], params.amplicon_info)

  if (params.QC_only == false) {
    
    DADA2_ANALYSIS(CUTADAPT.out[0], params.amplicon_info)

    if (params.refseq_fasta == null || params.masked_fasta == null) {
      if (params.genome == null) {
        error 1, "If reference sequences are not provided, a path to a genome must be provided to create reference sequences"
      }

      CREATE_REFERENCE_SEQUENCES(params.amplicon_info, params.genome, "${params.target}_refseq.fasta")
      DADA2_POSTPROC(DADA2_ANALYSIS.out[0], params.homopolymer_threshold, CREATE_REFERENCE_SEQUENCES.out[0], CREATE_REFERENCE_SEQUENCES.out[1], params.parallel)
    
      if (params.genome == null) {
        error 1, "If reference sequences are not provided, a path to a genome must be provided to create reference sequences"
      }

      // BUILD_GENOME_INDEX(params.genome)
      RESISTANCE_MARKERS(DADA2_POSTPROC.out[0], CREATE_REFERENCE_SEQUENCES.out[0], params.genome, codontable, params.resmarkers_amplicon)

    } else {
      DADA2_POSTPROC(DADA2_ANALYSIS.out[0], params.homopolymer_threshold, params.refseq_fasta, params.masked_fasta, params.parallel)

      if (params.genome == null) {
        error 1, "If reference sequences are not provided, a path to a genome must be provided to create reference sequences"
      }

      // BUILD_GENOME_INDEX(params.genome)
      RESISTANCE_MARKERS(DADA2_POSTPROC.out[0], CREATE_REFERENCE_SEQUENCES.out[0], params.genome, codontable, params.resmarkers_amplicon)

    }
  }
}

process CREATE_REFERENCE_SEQUENCES {

          input:
          path amplicon_info
          path genome
          val refseq_fasta

          output:
          path "*.fasta"
          path "*.mask"

          script:

          """
          Rscript ${params.scriptDIR}/create_refseq.R ${amplicon_info} ${genome} ${refseq_fasta}
          trf "${params.target}_refseq.fasta" 2 7 7 80 10 25 3 -h -m
          """
}


// cutadapt based quality filtering and QC
// Filter to only reads with primer detected and trim poly-g
process CUTADAPT {

        conda (params.enable_conda ? 'envs/cutadapt-env.yml' : null)

        input:
        tuple val(pair_id), path(reads)
        path amplicon_info
        val cutadapt_minlen
        val qualfilter

        output:
        file('trimmed_demuxed')
        file '*{AMPLICON,SAMPLE}summary.txt'

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
process QUALITY_CHECK {

        conda (params.enable_conda ? 'pandoc' : null)

        publishDir "${outDIR}",
               saveAs: { filename -> "${filename}"
               }

        input:
        tuple val(pair_id), path(reads) 
        path summfile
        path amplicon_info

        output:
        path '*coverage.txt'
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

        test -d quality_report || mkdir quality_report
        Rscript ${params.scriptDIR}/cutadapt_summaryplots.R amplicon_coverage.txt sample_coverage.txt ${amplicon_info} quality_report
        """
}

// Dada2

process DADA2_ANALYSIS {

        publishDir "${outDIR}",
               saveAs: { filename -> "${filename}"
               }

        input:
        path 'trimmed_demuxed'
        path amplicon_info

        output:
        path '*.RData'
        
        script:
        treat_no_overlap_differently = 'T'

        """
        Rscript ${params.scriptDIR}/dada_overlaps.R ${trimmed_demuxed} ${amplicon_info} $treat_no_overlap_differently dada2_output.RData
        """
}

// Dada2 Postprocessing
process DADA2_POSTPROC {

        publishDir "${outDIR}",
               saveAs: { filename -> "${filename}"
               }

        input:
        path rdatafile
        val homopolymer_threshold
        path refseq_fasta
        path masked_fasta
        val parallel

        output:
        path '*.{RDS,txt,csv}'
        
        script:
          """
          Rscript ${params.scriptDIR}/postdada_rearrange.R \
            --dada2-output ${rdatafile} \
            --homopolymer-threshold ${homopolymer_threshold} \
            --refseq-fasta ${refseq_fasta} \
            --masked-fasta ${masked_fasta} \
            --parallel ${parallel}
          """
}


// Resistance markers
process RESISTANCE_MARKERS {

        publishDir "${params.outDIR}",
               saveAs: {filename -> "${filename}"
               }

 input:
        path asvfile
        path refseq_fasta
        path genome
        path codontable
        path resmarkers_amplicon

 output:
        path '*.txt' 
        file('Mapping')
             
        script:
        """
       #!/usr/bin/env bash
       set -e
       bwa index ${genome}
       rdsfile2="\$(echo $asvfile | tr ' ' '\n' | grep -v RDS)"
       echo "\${rdsfile2}"
       mkdir -p Mapping
       cp ${refseq_fasta} Mapping
        awk 'BEGIN{FS=OFS="\\t";} {if(NR !=1) {print \$5,\$3}}' "\${rdsfile2}" | sort -u | awk 'BEGIN{FS=OFS="\\t"}{print">"\$1"\\n"\$2 >"Mapping/"\$1".fa"}'


       cd Mapping
        for bfile in *.fa; do allele=`echo \$bfile | cut -f 1-2 -d '.'`; echo \$allele; bwa mem -L 10000 ../${genome} \$bfile | samtools sort -o \$allele".bam" - ; samtools index \$allele".bam" ; done

        for cfile in *.bam; do allele=`echo \$cfile | cut -f 1-2 -d '.'`; echo \$allele; 
bcftools mpileup -d 2000 -f ../${refseq_fasta} \$cfile | bcftools query --format '%CHROM\\t%POS\\t%REF\\t%ALT\\n' > \$allele".mpileup.txt"; done
       cd ..
       Rscript ${params.scriptDIR}/resistance_marker_genotypes_bcftools_v4.R "\${rdsfile2}" ${codontable} ${resmarkers_amplicon} Mapping
       
        """
}

