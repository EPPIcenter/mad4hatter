# AmpSeq Workflow

## Setup

### Setting Parameters

Modify the nextflow.config file:

|Parameter|Description|
|---|---|
|readDIR|The folder that contains all the fastq files (*required*)|
|outDIR|The folder where you want the resulting data to be save (default 'results')|
|sequencer|The sequencer used to produce your data (default 'nextseq')|
|QC_only|Whether to only run QC related workflows or all workflows|
|refseq_fasta **or** genome|Path to reference sequences **or** path to genome (*one* is **required**)|

Additionally, the nextflow parameter `-profile` can be use to target the infrastructure you wish to run the pipeline on. The different profiles are listed below, including any setup that is required.

### Singularity

If using singularity, please run the command below to generate the singularity image.

```bash
sudo singularity build ampseq_worfklow.sif Singularity
```

And then include the `singularity` profile on the command line. 

*Note: you should also include executor you wish to run*

```bash
nextflow run main.nf --readDIR single --refseq_fasta v4_refseq.fasta --target v4 -profile sge,singularity
```

Below is an example using the genome parameter:

```bash
nextflow run main.nf --readDIR ~/Documents/MAD4HATTER_example_data/single -w ~/Documents/work --target v4 -profile sge,singularity --genome PlasmoDB-59_Pfalciparum3D7_Genome.fasta
```

### Docker

The pipeline can be easily run with docker and is the recommended way to run it when not using an HPC.

Follow the steps below to setup your docker image:

*Note: [docker](https://www.docker.com/) is a prerequisite.*

```bash
docker build -t aarandad/ampseq_worfklow .
```

And you're done! To run the pipeline, simply add `-profile docker`. 

```bash
nextflow run main.nf --readDIR single --target v4-profile -profile docker
```

### Conda

To use conda, you must first install either [conda](https://docs.conda.io/en/latest/) or [miniconda](https://docs.conda.io/en/latest/miniconda.html). Once installed, include the `conda` profile on the command line.

```bash
nextflow run main.nf --readDIR single --target v3 -profile conda
```

### Customizing for your institution

There is a file named `custom.config` in `conf/` that can be used to tailor processes to your environment. By default,
this file is used to tailor the pipeline for Wynton HPC at UCSF. This file may be altered to fit your institution's profile.

### Examples 

Potential ways to execute the pipeline:

```bash
# local executor
nextflow run main.nf --target v3

# with a profile (currently only supports sge)
nextflow run main.nf -profile sge --target v3

# run with singularity on an HPC with specified reference sequences
nextflow run main.nf --readDIR single -profile sge,singularity --refseq_fasta v4_refseq.fasta --target v4

# or locally with docker
nextflow run main.nf --readDIR ~/Documents/MAD4HATTER_example_data/single/ --target v4 -profile docker --refseq_fasta v4_refseq.fasta

# genomes can be provided in lieu of reference sequences, which will be generated with the amplicon table
nextflow run main.nf --readDIR ~/Documents/MAD4HATTER_example_data/single/ -w ~/work --target v4 -profile docker --genome PlasmoDB-59_Pfalciparum3D7_Genome.fasta
```

If you need to resume after some processes were successfully executed, add -resume at the end of it

