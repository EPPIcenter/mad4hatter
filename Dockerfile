FROM rocker/tidyverse:4.3.1

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    python3-pip \
    libssl-dev \
    libbz2-dev \
    libsdl1.2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    trf \
    bwa \
    bcftools \
    samtools \
    default-jre \
    pandoc \
    jq \
    parallel \
    fastp \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies
RUN pip3 install cutadapt==4.4 pandas==2.0.3 biopython==1.79

# R configuration
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'

# R packages
RUN Rscript -e 'install.packages("remotes")'
RUN Rscript -e 'remotes::install_cran(c("ggbeeswarm", "gridExtra", "rmarkdown", "foreach", "doMC", "argparse"), upgrade = "never")'

# Bioconductor packages
RUN Rscript -e 'if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager"); }; BiocManager::install(version = "3.18");'
RUN Rscript -e 'BiocManager::install("dada2", ask = FALSE)'
RUN Rscript -e 'BiocManager::install("muscle", ask = FALSE)'
RUN Rscript -e 'BiocManager::install("BSgenome", ask = FALSE)'

# Install additional tools needed for the analysis
RUN Rscript -e 'BiocManager::install(c("ShortRead", "GenomicRanges", "IRanges"), ask = FALSE)'

# Ensure the locale is set to UTF-8
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
