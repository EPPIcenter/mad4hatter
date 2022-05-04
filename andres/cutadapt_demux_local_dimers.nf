// Script parameters
fwd_primers = file( params.fwd_primers )
rev_primers = file( params.rev_primers )
amplicon_info = file( params.amplicon_info )

cutadapt_minlen = params.cutadapt_minlen
if ( params.sequencer == 'miseq' ) { qualfilter = 10 } else { qualfilter = 20 }


/*
Create 'read_pairs' channel that emits for each read pair a
tuple containing 3 elements: pair_id, R1, R2
*/
Channel
    .fromFilePairs( params.reads )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { read_pairs } 

//Split the Channel
read_pairs.into { read_pairs_cutadapt; read_pairs_filter_dimers }


// cutadapt based quality filtering and QC
// Filter to only reads with primer detected and trim poly-g
process cutadapt {
    
        publishDir "${params.outDIR}",
                saveAs: {filename ->
                        if (filename.endsWith("cutadapt_summary.txt")) "${pair_id}/logs/${filename}"
                        else "${pair_id}/${filename}"
                }

        input:
        tuple pair_id, file(reads) from read_pairs_cutadapt
        file fwd_primers
        file rev_primers
        val cutadapt_minlen
        val qualfilter
        
        output:
        file("trimmed_demuxed") 
        tuple SOMETHING, file(SOMETHING_R1), file(SOMETHING_R2) into filter_dimers


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
            -j 2 \
            --minimum-length ${cutadapt_minlen} \
            -o trimmed_demuxed/{name}_${pair_id}_trimmed_R1.fastq.gz \
            -p trimmed_demuxed/{name}_${pair_id}_trimmed_R2.fastq.gz \
            ${reads[0]} \
            ${reads[1]} 1> ${pair_id}.cutadapt_summary.txt
        """
}
    
process filter_dimers {
        publishDir "${params.outDIR}"

        input:
        tuple SOMETHING from filter_dimers.collect()
         
        output:
        file("filtered_out")
        file("filtered_in")

        script:
        """
        #!/usr/bin/env bash
        set -e

        mkdir filtered_out
        mkdir filtered_in

        echo ${trimmed_demuxed}
        cutadapt \
            --action=retain \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -j 2 \
            -o filtered_out/{name}_trimmed_filtered_R1.fastq.gz \
            -p filtered_out/{name}_trimmed_filtered_R2.fastq.gz \
            --untrimmed-output filtered_in/{name}_trimmed_filtered_R1.fastq.gz \
            --untrimmed-paired-output filtered_in/{name}_trimmed_filtered_R2.fastq.gz \
            ${pair_id}_trimmed_R1.fastq.gz \
            ${pair_id}_trimmed_R1.fastq.gz
        """

}

