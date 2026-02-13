# MAD4HATTER Amplicon Sequencing Pipeline

MAD4HATTER is a bioinformatics analysis pipeline used to process amplicon sequencing data. 

A quick start guide for running the complete pipeline is below and comprehensive documentation can be found [here](https://eppicenter.github.io/wonderland-docs/).

## Contents

- [Setup](#setup)
- [Setting Parameters](#setting-parameters)
    - [Mandatory Parameters](#mandatory-parameters)
    - [Optional Parameters](#optional-parameters)
        - [DADA parameters](#dada-parameters)
        - [Post processing parameters](#post-processing-parameters)
        - [Resmarker Module Parameters](#resmarker-module-parameters)
- [Runtime Profiles](#runtime-profiles)
    - [Apptainer](#apptainer)
    - [Docker](#docker)
    - [Conda](#conda)

## Setup

The mad4hatter pipeline uses [nextflow](https://www.nextflow.io/) and this will need to be installed prior to using the pipeline. Information about how to [install](https://www.nextflow.io/) and use the [command line tool](https://www.nextflow.io/docs/latest/cli.html) can be found on their [website](https://www.nextflow.io/). The tool is also available from other package managers such as [conda](https://anaconda.org/bioconda/nextflow) if you would like an alternative installation pathway. 

One of the useful features of Nextflow is that it caches your job history, so if for any reason your pipeline fails midway you can make changes to fix the failure and use the `-resume` flag to pick up where you left off. See more information [here](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

### Setting Parameters

To view the parameters and examples on the command line, run:

```bash 
nextflow run main.nf --help 
```

#### Mandatory Parameters

Below are the parameters that are essential for running the complete pipeline. For more details on running just a portion of the pipeline steps [see here](https://eppicenter.github.io/wonderland-docs/execution/).
 
|Parameter|Description|
|---|---|
|pools|The pools that were used for sequencing. [Options: D1,R1,R2 - check panel.config for more options]|
|readDIR|Path to folder containing fastq files|

Here is an example of running the complete workflow: 

```bash
nextflow run main.nf --readDIR /path/to/data --pools D1,R1,R2 -profile sge,apptainer
``` 

#### Optional Parameters

Below are parameters that are optional for running the pipeline.

|Parameter|Description|
|---|---|
| outDIR | The folder where you want the resulting data to be saved (default 'results') |
| workflow_name | Workflow option to be run [Options: complete (default), qc, postprocessing] |
|**Nextflow parameters**||
|profile|The infrastructure you wish to run the pipeline on. The different profiles are listed below under `Runtime Profiles`, including any setup that is required. **Please read that section for more details.**|
|config|Resource configurations for each process that will override any defaults set in the pipeline. It is recommended to use the provided `custom.config` file to make these resource modifications.|

Below is an example of how you may run the pipeline setting the above parameters. 

```bash
nextflow run main.nf --readDIR /path/to/data --outDIR /path/to/results --pools D1,R1,R2 -profile docker --workflow_name qc -config conf/custom.config 
``` 

##### DADA parameters

DADA2 infers amplicon sequences exactly and can be tuned depending on your needs. DADA2 is run in the DADA2 module of the pipeline (`DADA2_ANALYSIS`). Below are parameters that you can set to control your output. 

|Parameter|Description|
|---|---|
|omega_a|This controls the abundance threshold used to determine whether a sequence is overly abundant such that it is likely a true variant and not an error produced by DADA. (default `1e-120`)|
|dada2_pool|The method for information sharing across samples (default `pseudo`)|
|band_size|An alignment heuristic that controls whether an alignment will occur between sequences if the number of indels exceeds this threshold (default `16`)|
|maxEE|During filtering and trimming, reads that exceed the number of expected errors will be discarded (default `3`)|
|just_concatenate|Setting this to true will concatenate any DADA sequences that were unable to be merged. Reads that are concatenated will have 10 Ns separating the forward and reverse reads (i.e. `N`). Setting this to false will discard reads that did not have enough bases to merge. The minimum overlap required to merge forward and reverse reads is 12 bases. (default true)|
TODO: come back to this 

For more information about DADA2 and the parameters that can be set, please refer to their [documentation](https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf). 

Below is an example of how you may use the above parameters on the command line:

```bash
nextflow run main.nf --readDIR /path/to/data --outDIR /path/to/results -profile docker --pools D1,R1,R2 -config conf/custom.config --omega_a 1e-120 --band_size 16 --dada2_pool pseudo
```

##### Post processing parameters

By default the pipeline will use the `--pools` parameter and `panel.config` to find the paths to the reference sequences for each of the pools. This can be overridden by setting either the `refseq_fasta` **or** `genome` parameter as detailed below. 

Below are parameters that you can set to control the postprocessing module.

|Parameter|Description|
|---|---|
|refseq_fasta **or** genome|Path to targeted reference sequence **or** a specified genome that covers all targets. If neither are specified then a reference will be built from the fasta files under `panel_information` based on the pools supplied.|
|homopolymer_threshold|Homopolymers greater than this threshold will be masked (default `5`)|
|trf_min_score|Used by Tandem Repeat Finder. This will control the alignment score required to call a sequence a tandem repeat and mask it (default `25`)|
|trf_max_period|Used by Tandem Repeat Finder. This will limit the range of the pattern size of a tandem repeat to be masked (default `3`)|

Below is a continuation of the example above that shows how these parameters may be modified on the command line. Note that `--refseq_fasta` OR `--genome` OR no flag can be set to provide a reference. If no reference is provided (neither `--refseq_fasta` OR `--genome` are set) then the pipeline will build a targeted reference from the reference for each pool, stored under the resources directory. 

```bash
nextflow run main.nf --readDIR /path/to/data --outDIR /path/to/results -profile sge,apptainer --refseq_fasta /path/to/targeted_reference --pools D1,R1,R2 -config conf/custom.config --omega_a 1e-120 --band_size 16 --dada2_pool pseudo --trf_min_score 25 --trf_max_period 3
```

```bash
nextflow run main.nf --readDIR /path/to/data --outDIR /path/to/results -profile sge,apptainer --genome /path/to/Whole_Genome.fasta --pools D1,R1,R2 -config conf/custom.config --omega_a 1e-120 --band_size 16 --dada2_pool pseudo --trf_min_score 25 --trf_max_period 3
```
##### Resmarker Module Parameters 

By default the pipeline will check if any of the markers in the principal list (`panel_information/principal_resistance_marker_info_table.tsv`) are covered by any of the loci in the panel. If markers are covered then the resmarker module will run. This can be overridden and a customized table can be supplied using the following parameter.

|Parameter|Description|
|---|---|
|resmarker_info|Path to table containing resmarker information|

## Runtime Profiles

Runtime profiles will provide all dependencies and setup needed for different computing environments. As an example, if you are using a cluster, grid or HPC environment, `apptainer` would be an appropriate profile as it supplies an image with all dependencies ready. If you are using a local computer, `docker` would be more appropriate. You can also choose to install the dependencies independently and run the pipeline that way if you choose to, but it is not recommended. 

Continuing with our example above, the below could be used to run the pipeline using the SGE scheduler.

```bash
nextflow run main.nf --readDIR /path/to/data -profile sge,apptainer --pools D1,R1,R2 
```

### Apptainer

*Note: [apptainer](https://github.com/apptainer/apptainer/releases) is a prerequisite.*

Apptainer should be used if you are using a computing cluster or grid. All dependencies needed to run the pipeline are contained within the apptainer image. The image can be created by pulling the docker image from dockerhub, which will create a `mad4hatter_latest.sif` image in your working directory.

```bash
apptainer pull docker://eppicenter/mad4hatter:latest
```

As an alternative, you can build the image yourself by running the command below.

```bash
apptainer build mad4hatter_latest.sif Apptainer
```

Once you have the image, you must include the `apptainer` profile on the command line in order for it to be used.

*Note: you should also include the job scheduler you will be using. In this case, `sge` is the job scheduler that will be used. Contact your system administrator if you are unsure about this setting.*

```bash
nextflow run main.nf --readDIR single --pools D1,R1,R2 -profile sge,apptainer
```

### Docker

*Note: [docker](https://www.docker.com/) is a prerequisite.*

The pipeline can be easily run with docker and is the recommended way to run it when not using an HPC.

The EPPIcenter has a repository for images, and the docker image for the pipeline will be automatically pulled in the background when first running the pipeline. The image will then be stored locally on your machine and reused. 

To run the pipeline with Docker, simply add `-profile docker` to your command. 

```bash
nextflow run main.nf --readDIR /path/to/data -profile docker --pools D1,R1,R2
```

Alternatively, you can build the docker image on your machine using the Dockerfile recipe, although this is not the recommended way to set up the docker image.

If you would like to build the docker image yourself, you may run the command below:

```bash
docker build -t eppicenter/mad4hatter:latest .
```

### Conda

To use conda, you must first install either [conda](https://docs.conda.io/en/latest/) or [miniconda](https://docs.conda.io/en/latest/miniconda.html). Once installed, include the `conda` profile on the command line.

```bash
nextflow run main.nf --readDIR /path/to/data -profile conda --pools D1,R1,R2
```