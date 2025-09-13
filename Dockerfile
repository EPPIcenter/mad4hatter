FROM rocker/tidyverse:4.3.1

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3-pip libssl-dev libbz2-dev libsdl1.2-dev liblzma-dev libcurl4-openssl-dev zlib1g-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev trf bwa bcftools samtools default-jre pandoc jq parallel && rm -rf /var/lib/apt/lists/*

# Python dependencies
RUN pip install numpy==1.25.2 cutadapt==4.4 pandas==2.0.3 bio==1.5.9 

# R configuration
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'

# R packages
RUN Rscript -e 'install.packages("remotes", version = "2.4.2")'
RUN Rscript -e 'remotes::install_cran(c("ggbeeswarm", "gridExtra", "rmarkdown", "foreach", "doMC", "argparse"), upgrade = "never", version = c("0.6.1", "2.3", "2.17", "1.5.2", "1.3.8", "2.1.6"))'

# Bioconductor packages
RUN Rscript -e 'if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager"); }; BiocManager::install(version = "3.17", ask = FALSE);'
RUN Rscript -e 'BiocManager::install("dada2", version = "3.17", ask = FALSE)'
RUN Rscript -e 'BiocManager::install("muscle", version = "3.17", ask = FALSE)'
RUN Rscript -e 'BiocManager::install("BSgenome", version = "3.17", ask = FALSE)'
