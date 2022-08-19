# AmpSeq Workflow

## Setup

### Environment

```bash
conda env create -f environment.yaml
```

### R Dependencies

```R
install.packages(c('ggbeeswarm', 'tidyverse', 'gridExtra', 'rmarkdown', 'knitr'))

remotes::install_github("m-murphy/moire")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")
```

## Running the Pipeline

### Setting Parameters

Modify the nextflow.config file:

|Parameter|Description|
|---|---|
|readDIR|The folder that contains all the fastq files (*required*)|
|outDIR|The folder where you want the resulting data to be save (default 'results')|
|sequencer|May leave as miseq for now. Not functional at the moment (default 'nextseq')|
|QC_only|Whether to only run QC related workflows or all workflows (default 'F')|

### Customizing for your institution

There is a file named `custom.config` in `conf/` that can be used to tailor processes to your environment. By default,
this file is used to tailor the pipeline for Wynton HPC at UCSF. This file may be altered to fit your institution's profile.

### Executing the Pipeline

Potential ways to execute the pipeline:

```bash
# local executor
nextflow run main.nf

# with a profile (currently only supports sge)
nextflow run main.nf -profile sge
```

If you need to resume after some processes were successfully executed, add -resume at the end of it

