# AmpSeq Workflow

## Setup

### Environment

```bash
conda env create -f environment.yaml
```

### R Dependencies

```R
install.packages(c('dada', 'ggbeeswarm', 'tidyverse', 'gridExtra', 'rmarkdown', 'knitr'))
```

## Running the Pipeline

### Setting Parameters

To run copy main.nf (or main_wynton.nf) and corresponding config file to the folder where you have the data (you may also just call the code from that folder)
Modify in the config file:

|Parameter|Description|
|---|---|
|readDIR|The folder that contains all the fastq files|
|outDIR|The folder where you want the resulting data to be save|
|sequencer|May leave as miseq for now. Not functional at the moment|
|primerDIR|Path to the fastas folder|
|fwd_ and rev_ primers|Change to v3 or v4 depending on the version of the panel you're using|
|amplicon_info|Path to the info.tsv file. Use v3 or v4 as needed|
|scriptDIR|Path to the folder that has the R code.|
|process.conda|Path to your cutadapt environment|

### Executing the Script

Commands to run:

```bash
nextflow run main.nf -c nexflow.config -profile conda 
````

If you need to resume after some processes were successfully executed, add -resume at the end of it

