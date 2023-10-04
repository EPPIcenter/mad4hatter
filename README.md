# MAD4HATTER Amplicon Sequencing Pipeline

## Setup

The mad4hatter pipeline uses [nextflow](https://www.nextflow.io/) and will need to be installed prior to using the pipeline. Information about how to [install](https://www.nextflow.io/) and use the [command line tool](https://www.nextflow.io/docs/latest/cli.html) can be found on their [website](https://www.nextflow.io/). The tool is also available from other package managers such as [conda](https://anaconda.org/bioconda/nextflow) if you would like an alternative installation pathway. 

### Setting Parameters

Below are parameters that are essential to running the pipeline. You will need to specify where your `fastq` files will be read from, what sequencer was used to generate the data, and what amplicon panel was used. You may optionally include a destination folder where you would like final pipeline output to go, or choose to only run the QC portion of the pipeline. 

|Parameter|Description|
|---|---|
|readDIR|The folder that contains all the fastq files (*required*)|
|outDIR|The folder where you want the resulting data to be save (default 'results')|
|sequencer|The sequencer used to produce your data (default 'nextseq')|
|QC_only|Whether to only run QC related workflows or all workflows|
|target|The amplicon panel that was used to create your sequencing data (*required*)|
|profile|The infrastructure you wish to run the pipeline on. The different profiles are listed below under `Runtime Profiles`, including any setup that is required. **Please read that section for more details.**|
|config|Resource configurations for each process that will override any defaults set in the pipeline. It is recommend to use the provided `custom.config` file to make these resource modifications.|

Below is an example of how you may run the pipeline setting the above parameters. 

```bash
nextflow run main.nf --readDIR /wynton/scratch/data --outDIR /wynton/scratch/results --target v4 -profile sge,apptainer --sequencer miseq --QC_only -config conf/custom.config 
``` 

### DADA parameters

DADA2 infers amplicon sequences exactly and can be tuned depending on your needs. DADA2 is run in the DADA2 module of the pipeline (`DADA2_ANALYSIS`). Below are parameters that you can set to control your output. 

|Parameter|Description|
|---|---|
|omega_a|This controls the abundance threshold used to determine whether a sequence is overly abundant such that it is likely a true variant and not an error produced by DADA. (default `1e-120`)|
|pool|The method for information sharing across samples (default `pseudo`)|
|band_size|An alignment heursitic that controls whether an alignment will occur between sequences if the number of indels exceed this threshold (default `16`)|
|maxEE|During filtering and trimming, reads that exceed the number of expected errors will be discarded (default `2`)|
|concat_non_overlaps|Setting this to true will concatenate any DADA sequences that were unable to be merged. Reads that are concatenated will have 10 Ns separating the forward and reverse read (ie. `NNNNNNNNNN`) Setting this to false will discard reads that did not have enough bases to merge. The minimum overlap required to merge forward and reverse reads is 12 bases.|

For more information about DADA2 and the parameters that can be set, please refer to their [documentation](https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf). 

Below is an example of how you may use the above parameters on the command line:

```bash
nextflow run main.nf --readDIR /wynton/scratch/data --outDIR /wynton/scratch/results -profile sge,apptainer --target v4 -config conf/custom.config --omega_a 1e-120 --band_size 16 --pool pseudo
```

### Post processing parameters

There are a number of steps in the postprocessing module (`DADA2_POSTPROC`) to reduce false positives in your final allele table. In particular, tandem repeats and hompolymer regions are masked due to known problems that they cause in Illumina instruments. Masked regions will appear as `N`'s in your sequences such as the example below.

Here is a sequence from DADA:
```
TATATATATATATATATATATATATATATATATATATATATATATGTATGTATGTTGATTAATTTGTTTATATATTTATATTTATTTCTTATGACCTTTTTAGGAACGACACCGAAGCTTTAATTTACAATTTTTTGCTATATCCATGTTAGATGCCTGTTCAGTCATTTTGGCCTTCATAGGTCT
```

And here is it's masked counterpart. Notice how the homopolymers and tandem repeats are masked by `N`s.
```
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGTATGTATGTTGATTAATTTGTTTATATATTTATATTTATTTCTTATGACCTTTTTAGGAACGACACCGAAGCTTTAATTTACANNNNNNNNCTATATCCATGTTAGATGCCTGTTCAGTCATTTTGGCCTTCATAGGTCT
```

The masking is accomplished using [Tandem Repeat Finder](https://github.com/Benson-Genomics-Lab/TRF#trf-definitions). Please refer to their documentation for additional information.

Below are parameters that you can set to control the postprocessing module.

|Parameter|Description|
|---|---|
|refseq_fasta **or** genome|Path to reference sequences **or** a specified genome target that is registered in the pipeline. Additionally you can specify a path to genome (*one* is **required**). This is required in order to map DADA sequences to their respective amplicons and to identify off target sequences|
|add_mask|Whether to add a mask or not to the final sequences (default `true`)|
|trf_min_score|Used by Tandem Repeat Finder. This will control the alignment score required to call a sequence a tandem repeat (default `25`)|
|trf_max_period|Used by Tandem Repeat Finder. This will limit the range of the pattern size of a tandem repeat (default `3`)|

Below is a continuation of the example above that shows how these parameters may be modified on the command line. Notice the `genome` parameter will accept an identifier or a path to a fasta. Details about the fasta associated with an identifier can be found in `conf/base.config`.

```bash
nextflow run main.nf --readDIR /wynton/scratch/data --outDIR /wynton/scratch/results -profile sge,apptainer --genome v1 --target v4 -config conf/custom.config --omega_a 1e-120 --band_size 16 --pool pseudo --trf_min_score 25 --trf_max_period 3
```

```bash
nextflow run main.nf --readDIR /wynton/scratch/data --outDIR /wynton/scratch/results -profile sge,apptainer --genome /wynton/share/PlasmoDB-59_Pfalciparum3D7_Genome.fasta --target v4 -config conf/custom.config --omega_a 1e-120 --band_size 16 --pool pseudo --trf_min_score 25 --trf_max_period 3
```

## Runtime Profiles

Runtime profiles will provide all dependencies and setup needed for different computing environments. As an example, if you are using a cluster, grid or HPC environment, `apptainer` would be an appropriate profile as it supplies an image with all dependencies ready. If you are using a local computer, `docker` would be more appropriate. You can also choose to install the dependencies independently and run the pipeline that way if you choose to, but it is not recommended. 

Currently, [Sun Grid Engine (SGE)](http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html) is the only supported cluster environment.

Continuing with our example above, the below could be used to run the pipeline using the SGE scheduler.

```bash
nextflow run main.nf --readDIR /wynton/scratch/data --outDIR /wynton/scratch/results -profile sge,apptainer --genome /wynton/share/PlasmoDB-59_Pfalciparum3D7_Genome.fasta --target v4 -config conf/custom.config --omega_a 1e-120 --band_size 16 --pool pseudo --trf_min_score 25 --trf_max_period 3
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
nextflow run main.nf --readDIR single --refseq_fasta v4_refseq.fasta --target v4 -profile sge,apptainer -c conf/custom.config
```

### Docker

*Note: [docker](https://www.docker.com/) is a prerequisite.*

The pipeline can be easily run with docker and is the recommended way to run it when not using an HPC.

The EPPIcenter has a repository for images, and the docker image for the pipeline will be automatically pulled in the background when first running the pipeline. The image will then be stored locally on your machine and reused. 

To run the  with docker, simply add `-profile docker` in your command. 

```bash
nextflow run main.nf --readDIR /wynton/scratch/data --outDIR /wynton/scratch/results -profile docker --genome /wynton/share/PlasmoDB-59_Pfalciparum3D7_Genome.fasta --target v4 -config conf/custom.config
```

Alternatively, you can build the docker image on your machine using the Dockerfile recipe, although this is not the recommended way to set up the docker image.

If you would like to build the docker image yourself, you may run the command below:

```bash
docker build -t eppicenter/mad4hatter:latest .
```

### Conda

To use conda, you must first install either [conda](https://docs.conda.io/en/latest/) or [miniconda](https://docs.conda.io/en/latest/miniconda.html). Once installed, include the `conda` profile on the command line.

```bash
nextflow run main.nf --readDIR /wynton/scratch/data --outDIR /wynton/scratch/results -profile conda --genome /wynton/share/PlasmoDB-59_Pfalciparum3D7_Genome.fasta --target v4 -config conf/custom.config
```
