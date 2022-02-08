// Script parameters
fwd_primers = file( params.fwd_primers )
rev_primers = file( params.rev_primers )

params.dev = false
params.number_of_inputs = 1

/*
Create 'read_pairs' channel that emits for each read pair a
tuple containing 3 elements: pair_id, R1, R2
*/
Channel
    .fromFilePairs( params.reads )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .take( params.dev ? params.number_of_inputs : -1 )
        .set { read_pairs }

// cutadapt based quality filtering and QC
// Filter to only reads with primer detected and trim poly-g
process cutadapt {
    
        publishDir "${params.outdir}",
                saveAs: {filename ->
                        if (filename.indexOf(".log") > 0) "${pair_id}/logs/${filename}"
                        else if (filename.indexOf(".html") > 0) "${pair_id}/htmls/${filename}"
                        else "${pair_id}/${filename}"
                }

        input:
        set pair_id, file(reads) from read_pairs
        file fwd_primers
        file rev_primers

        output:
        file("trimmed_demuxed")

        conda 'bioconda::cutadapt=3.5'

        time '120m'
        cpus 8
        penv 'smp'
        memory '8 GB'

        script:
        """
        #!/usr/bin/env bash
        set -e
        ulimit -s unlimited

        mkdir trimmed_demuxed

        cutadapt \
            --action=trim \
            --discard-untrimmed \
            -g file:${fwd_primers} \
            -G file:${rev_primers} \
            --pair-adapters \
            -e 0 \
            --no-indels \
            -o trimmed_demuxed/{name}_${pair_id}_trimmed_R1.fastq.gz \
            -p trimmed_demuxed/{name}_${pair_id}_trimmed_R2.fastq.gz \
            ${reads[0]} \
            ${reads[1]} 
            
        """
}
