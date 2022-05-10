# ampseq_workflow

Install:
R (will need libraries: dada, ggbeeswarm, tidyverse)
cutadapt
nextflow

To run copy main.nf (or main_wynton.nf) and corresponding config file to the folder where you have the data (you may also just call the code from that folder)
Modify in the config file:

readDIR: the folder that contains all the fastq files
outDIR: the folder where you want the resulting data to be saved
sequencer: may leave as miseq for now. Not functional at the moment
primerDIR: path to the fastas folder
fwd_ and rev_ primers: change to v3 or v4 depending on the version of the panel you're using
amplicon_info: path to the info.tsv file. Use v3 or v4 as needed 
scriptDIR: path to the folder that has the R code. 
process.conda: path to your cutadapt environment


Commands to run:

nextflow run main.nf -c nexflow.config -profile conda 

If you need to resume after some processes were successfully executed, add -resume at the end of it

