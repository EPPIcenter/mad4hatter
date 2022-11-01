# AmpSeq Workflow

## Setup

### R Dependencies

```R
install.packages(c('ggbeeswarm', 'tidyverse', 'gridExtra', 'rmarkdown', 'knitr'))

remotes::install_github("EPPIcenter/moire")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")
```

## Running the Pipeline

### Docker

The pipeline can be easily run with docker and is the recommended way to run it when not using an HPC.

Follow the steps below to setup your docker image:

*Note: [docker] (https://www.docker.com/) is a prerequisite.*

```bash
docker build -t aarandad/ampseq_worfklow .
```

And you're done! To run the pipeline, simply add `-profile docker`. 

```bash
nextflow run main.nf --readDIR </path/to/fastqs> -profile docker --target v4
```

### Setting Parameters

Modify the nextflow.config file:

|Parameter|Description|
|---|---|
|readDIR|The folder that contains all the fastq files (*required*)|
|outDIR|The folder where you want the resulting data to be save (default 'results')|
|sequencer|The sequencer used to produce your data (default 'nextseq')|
|QC_only|Whether to only run QC related workflows or all workflows (default 'F')|

### Customizing for your institution

There is a file named `custom.config` in `conf/` that can be used to tailor processes to your environment. By default,
this file is used to tailor the pipeline for Wynton HPC at UCSF. This file may be altered to fit your institution's profile.

### Executing the Pipeline

Potential ways to execute the pipeline:

```bash
# local executor
nextflow run main.nf --target v3

# with a profile (currently only supports sge)
nextflow run main.nf -profile sge --target v3
```

If you need to resume after some processes were successfully executed, add -resume at the end of it

