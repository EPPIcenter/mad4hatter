# Getting Started

This guide will help you get started, even if you're new to bioinformatics pipelines.

## Prerequisites

Before you begin, you'll need:

1. **Nextflow** - The workflow management system that runs the pipeline
    - Installation instructions: [Nextflow website](https://www.nextflow.io/)
    - Alternative: Install via [conda](https://anaconda.org/bioconda/nextflow): `conda install -c bioconda nextflow`

2. **Java 11 or higher** - Required by Nextflow
    - Check if you have Java: `java -version`
    - If not installed, download from [Oracle](https://www.oracle.com/java/technologies/downloads/) or use your system's package manager

3. **One of the following runtime environments** (choose based on your setup):
    - **Docker** - For local computers (recommended for beginners)
    - **Apptainer/Singularity** - For HPC clusters
    - **Conda** - Alternative dependency management

## Quick Start Guide

### Step 1: Choose Your Environment

- **Are you running on:**
    - **Your local computer?** → Use [Docker](#docker)
    - **A computing cluster/HPC?** → Use [Apptainer](#apptainer)
    - **Prefer conda?** → Use [Conda](#conda-alternative-option)

### Step 2: Prepare Your Data

Make sure your sequencing data is organized:
- Forward reads (R1) and reverse reads (R2) in FASTQ format
- Files should be in a single folder
- Example: `/path/to/data/sample1_R1.fastq.gz` and `/path/to/data/sample1_R2.fastq.gz`

### Step 3: Run the Pipeline

Follow the instructions for your chosen environment below.

---

## Docker

**Best for:** Running on your own computer (Mac, or Linux). Docker automatically handles all software dependencies. The pipeline will download the required Docker image automatically on first run - no manual setup needed!


### Prerequisites
- [Docker Desktop](https://www.docker.com/products/docker-desktop/) installed and running

### Basic Command

```bash
nextflow run main.nf \
  --readDIR tests/example_data/example_fastq \
  --pools D1,R1,R2 \
  -profile docker
```

!!! tip "First Time Running?"
    The first time you run with Docker, it will download the pipeline image (this may take a few minutes). Subsequent runs will be much faster!

---

## Apptainer

**Best for:** High-performance computing clusters, grid computing, or shared computing resources

### Prerequisites
- [Apptainer](https://github.com/apptainer/apptainer/releases) installed on your cluster
- Access to a cluster with a job scheduler (currently supports SGE or slurm)

### Step 1: Build the Apptainer Image

First, pull the container image onto your cluster:

```bash
apptainer pull docker://eppicenter/mad4hatter:latest
```

Alternatively, you can build the container image from scratch on your cluster:

```bash
apptainer build mad4hatter.sif Apptainer
```

This creates a file called `mad4hatter.sif` that contains all the software needed.

!!! note "One-Time Setup"
    You typically only need to build/pull the image once. After that, you can reuse the `mad4hatter.sif` file.

### Step 2: Run the Pipeline

```bash
nextflow run main.nf \
  --readDIR tests/example_data/example_fastq \
  --pools D1,R1,R2 \
  -profile sge,apptainer \
  -c conf/custom.config
```

**Important parameters:**
- `-profile sge,apptainer` - Use a job scheduler (e.g., SGE or slurm) with Apptainer
- `-c conf/custom.config` - Configuration file for resource allocation

!!! note "Job Scheduler"
    Replace `sge` with your cluster's job scheduler if different. Contact your system administrator if you're unsure which scheduler your cluster uses.

---

## Conda (Alternative Option)

**Best for:** Users who prefer conda for package management. 

!!! tip "Using conda"
    This is not a recommended way to run the pipeline and limitted support will be available for running using conda.

### Prerequisites
- [Conda](https://docs.conda.io/en/latest/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed

### Basic Command

```bash
nextflow run main.nf \
  --readDIR tests/example_data/example_fastq \
  --pools D1,R1,R2 \
  -profile conda
```

!!! note "Conda Environments"
    Conda will automatically create and manage the required software environments. This may take longer on the first run as it installs dependencies.

---

## No-Code Option: Terra Platform

**Don't want to use the command line?** The pipeline is also available on Terra, a cloud-based platform with a graphical interface.

- Access the workspace: [Terra MAD4HATTER Workspace](https://app.terra.bio/#workspaces/gates-malaria/Mad4Hatter)
- No installation required!

---

## Next Steps

1. **Review the outputs** - See [Pipeline Outputs](pipeline-outputs.md) for details
2. **Understand pipeline usage** - Check the [Running the Pipeline](execution.md) page for more advanced ways of running the pipeline. 

## Getting Help

- **Command-line help**: See all available parameters and options by running:
  ```bash
  nextflow run main.nf --help
  ```
  This will display the complete help message with all pipeline parameters and their descriptions.

- **Documentation**: Browse the other pages in this documentation
- **Getting in touch**: Report bugs, feature requests, and questions as an issue on the [GitHub repository](https://github.com/EPPIcenter/mad4hatter/issues). Alternatively, reach out to the EPPIcenter team (kathryn.murie@ucsf.edu).
