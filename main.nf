// Script parameters
fwd_primers = file( params.fwd_primers )
rev_primers = file( params.rev_primers )
amplicon_info = file( params.amplicon_info )

cutadapt_minlen = params.cutadapt_minlen
if ( params.sequencer == 'miseq' ) { qualfilter = '--trim-n -q 10' } else { qualfilter = '--nextseq-trim=20' }

//params.dev = false
//params.number_of_inputs = 1

/*
Create 'read_pairs' channel that emits for each read pair a
tuple containing 3 elements: pair_id, R1, R2
*/
Channel
        .fromFilePairs( params.reads )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
//    .take( params.dev ? params.number_of_inputs : -1 )
        .set { read_pairs }

//Split the Channel
read_pairs.into { read_pairs_cutadapt; read_pairs_qualitycheck }

// cutadapt based quality filtering and QC
// Filter to only reads with primer detected and trim poly-g
process cutadapt {
    
        publishDir "${params.outDIR}",
                saveAs: {filename ->
                        if (filename.endsWith("cutadapt_summary.txt")) "${pair_id}/logs/${filename}"
//                       else if (filename.endsWith(".html")) "${pair_id}/htmls/${filename}"
                        else "${pair_id}/${filename}"
                }

        input:
        tuple pair_id, file(reads) from read_pairs_cutadapt
        file fwd_primers
        file rev_primers
        val cutadapt_minlen
        val qualfilter
        
        output:
        file("trimmed_demuxed") into dada2
        file '*{AMPLICON,SAMPLE}summary.txt' into cutadapt_report

 //       conda 'cutadaptenv'
 //         conda '/home/isglobal.lan/ddatta/miniconda3/envs/cutadaptenv'
  





        script:
        """
        #!/usr/bin/env bash
       set -e

     
        mkdir trimmed_demuxed
        mkdir trimmed_demuxed_unknown
        mkdir filtered_out
        mkdir filtered_in
        
        cutadapt \
            --action=trim \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -j 2 \
            -e 0 \
            --no-indels \
            --minimum-length ${cutadapt_minlen} \
            -o filtered_out/${pair_id}_dimer_R1.fastq.gz \
            -p filtered_out/${pair_id}_dimer_R2.fastq.gz \
            --untrimmed-output filtered_in/${pair_id}_filtered_R1.fastq.gz \
            --untrimmed-paired-output filtered_in/${pair_id}_filtered_R2.fastq.gz \
            ${reads[0]} \
            ${reads[1]} 1> ${pair_id}.cutadapt1_summary.txt
        


           grep -E "Total read pairs processed:" ${pair_id}.cutadapt1_summary.txt | paste -s -d";\n"  | sed 's/Total read pairs processed://' | sed 's/ //g' | sed 's/,//g' | awk -v var1=${pair_id} 'BEGIN{FS=" ";OFS="\t"}{print var1, "Input",\$1};' > ${pair_id}.SAMPLEsummary.txt
           
          grep -E "Pairs discarded as untrimmed:" ${pair_id}.cutadapt1_summary.txt | paste -s -d";\n"  | sed 's/Pairs discarded as untrimmed://' | sed -r 's/[(].*//' | sed 's/ //g' | sed 's/,//g' |awk -v var1=${pair_id} 'BEGIN{FS=" ";OFS="\t"}{print var1, "Total filtered",\$1};' >> ${pair_id}.SAMPLEsummary.txt
           


        cutadapt \
            --action=trim \
            -g file:${fwd_primers} \
            -G file:${rev_primers} \
            --pair-adapters \
            -e 0 \
            --no-indels \
            -j 2 \
            ${qualfilter} \
            --minimum-length ${cutadapt_minlen} \
            -o trimmed_demuxed/{name}_${pair_id}_trimmed_R1.fastq.gz \
            -p trimmed_demuxed/{name}_${pair_id}_trimmed_R2.fastq.gz \
            --untrimmed-output trimmed_demuxed_unknown/${pair_id}_unknown_R1.fastq.gz \
            --untrimmed-paired-output trimmed_demuxed_unknown/${pair_id}_unknown_R2.fastq.gz \
            filtered_in/${pair_id}_filtered_R1.fastq.gz \
            filtered_in/${pair_id}_filtered_R2.fastq.gz 1> ${pair_id}.cutadapt2_summary.txt
           
           grep -E "Pairs written" ${pair_id}.cutadapt2_summary.txt | cut -f 2- -d ":" | sed -r 's/[(].*//' | sed 's/ //g' | sed 's/,//g' | awk -v var1=${pair_id} 'BEGIN{FS=" ";OFS="\t"}{print var1, "Final filtered",\$1};' >> ${pair_id}.SAMPLEsummary.txt
           
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

        publishDir "${params.outDIR}", 
               saveAs: {filename -> "${filename}"
               }

        input:
        tuple pair_id, file(reads) from read_pairs_qualitycheck
        file summfile from cutadapt_report.collect().ifEmpty([])
        file amplicon_info

        output:
        file '*.txt' into qualitycheck_report

        script:
        """
        #!/usr/bin/env bash
       set -e 
       
       echo $summfile | tr ' ' '\n' | grep 'filt.AMPLICONsummary.txt' | tr '\n' ' ' | xargs cat > ALL_filt.AMPLICONsummary.txt
       echo $summfile | tr ' ' '\n' | grep 'trim.AMPLICONsummary.txt' | tr '\n' ' ' | xargs cat > ALL_trim.AMPLICONsummary.txt
       echo $summfile | tr ' ' '\n' | grep 'SAMPLEsummary.txt' | tr '\n' ' ' | xargs cat > ALL_SAMPLEsummary.txt
       
      if [ -s ALL_filt.AMPLICONsummary.txt ]; then
      awk 'NR == FNR { key[\$1,\$2] = \$3; next } { \$3 = ((\$1,\$2) in key) ? key[\$1,\$2] : 0 };1' OFS="\t"  ALL_filt.AMPLICONsummary.txt ALL_trim.AMPLICONsummary.txt > ALL_filt.final.AMPLICONsummary.txt
      else
      awk 'BEGIN{FS=OFS="\t"} { print \$1,\$2,0; }'  ALL_trim.AMPLICONsummary.txt > ALL_filt.final.AMPLICONsummary.txt
      fi
      
       Rscript ${params.scriptDIR}/cutadapt_summaryplots.R ALL_filt.final.AMPLICONsummary.txt ALL_SAMPLEsummary.txt ${amplicon_info} ${params.outDIR}  
        """
}

// Dada2 
/*
process dada2_analysis {

        publishDir "${params.outDIR}",
               saveAs: {filename -> "${filename}"
               }

 input:
       file 'trimmed_demuxed' from dada2.collect().ifEmpty([])
       file amplicon_info
         
        output:
        file '*.RData' into dada2_summary

        script:
       treat_no_overlap_differently = "T"
        
        """
      Rscript ${params.scriptDIR}/dada_overlaps.R ${trimmed_demuxed} ${amplicon_info} $treat_no_overlap_differently dada2_overlaps.RData
        """
}

// Dada2 Postprocessing
process dada2_postproc {

        publishDir "${params.outDIR}",
               saveAs: {filename -> "${filename}"
               }

 input:
        file rdatafile from dada2_summary

        output:
        file '*.{RDS,txt}' into dada2_proc

        script:

        """
      Rscript ${params.scriptDIR}/postdada_rearrange.R $rdatafile
        """
}
*/
