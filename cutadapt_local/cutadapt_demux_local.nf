// Script parameters
fwd_primers = file( params.fwd_primers )
rev_primers = file( params.rev_primers )

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
    
        publishDir "${params.outdir}",
                saveAs: {filename ->
                        if (filename.endsWith("cutadapt_summary.txt")) "${pair_id}/logs/${filename}"
                        else "${pair_id}/${filename}"
                }

        input:
        tuple pair_id, file(reads) from read_pairs_cutadapt
        file fwd_primers
        file rev_primers

        output:
        file("trimmed_demuxed") 
        file '*qctempsummary.txt' into cutadapt_report

        conda '/wynton/home/eppicenter/aarandad/miniconda3/envs/cutadaptenv/'

        script:
        """
        #!/usr/bin/env bash
        set -e

        mkdir trimmed_demuxed

        cutadapt \
            --action=trim \
            -g file:${fwd_primers} \
            -G file:${rev_primers} \
            --pair-adapters \
            -e 0 \
            --no-indels \
	    --nextseq-trim=20 \
            -j 2 \
            --minimum-length  116\
            -o trimmed_demuxed/{name}_${pair_id}_trimmed_R1.fastq.gz \
            -p trimmed_demuxed/{name}_${pair_id}_trimmed_R2.fastq.gz \
            ${reads[0]} \
            ${reads[1]} 1> ${pair_id}.cutadapt_summary.txt
           
        
        grep -E "Read 1 with adapter:" ${pair_id}.cutadapt_summary.txt | paste -s -d";\n"  | sed 's/Read 1 with adapter://' | sed -r 's/[(].*//' | sed 's/ //g' | sed 's/,//g' \
           | awk -v var1=${pair_id} 'BEGIN{FS=" ";OFS="\t"}{print var1, "Total trimmed",\$1};' >> ${pair_id}.qctempsummary.txt


        grep -E "Total read pairs processed:" ${pair_id}.cutadapt_summary.txt | paste -s -d";\n"  | sed 's/Total read pairs processed://' | sed 's/ //g' | sed 's/,//g' \
           | awk -v var1=${pair_id} 'BEGIN{FS=" ";OFS="\t"}{print var1, "Input",\$1};' >> ${pair_id}.qctempsummary.txt

        """
}

process qualitycheck {

        publishDir "${params.outdir}", 
               saveAs: {filename -> "${filename}"
               }

        input:
        tuple pair_id, file(reads) from read_pairs_qualitycheck
        file summfile from cutadapt_report.collect().ifEmpty([])

        output:
        file '*.txt' into qualitycheck_report

        script:
        """
        #!/usr/bin/env bash
       set -e 
       cat $summfile > ALL_qctempsummary.txt
       """
}
