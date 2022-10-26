FROM bioconductor/bioconductor_docker:RELEASE_3_14
FROM rocker/tidyverse:4.2.1

RUN apt-get update && apt-get install -y build-essential python3-pip libbz2-dev && rm -rf /var/lib/apt/lists/*
RUN pip install cutadapt

RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'

# Update apt-get
RUN Rscript -e 'install.packages("remotes", version = "2.4.2")'
RUN Rscript -e 'remotes::install_cran("ggbeeswarm", upgrade="never", version="0.6.1")'
RUN Rscript -e 'remotes::install_cran("gridExtra",upgrade="never", version = "2.3")'
RUN Rscript -e 'remotes::install_cran("rmarkdown",upgrade="never", version = "2.17")'
RUN Rscript -e 'remotes::install_cran("gridExtra",upgrade="never", version = "2.3")'
RUN Rscript -e 'remotes::install_cran("foreach",upgrade="never", version = "1.5.2")'
RUN Rscript -e 'remotes::install_cran("doMC",upgrade="never", version = "1.3.8")'

RUN Rscript -e 'BiocManager::install("dada2", version = "3.15", ask = FALSE)'
RUN Rscript -e 'BiocManager::install("muscle", ask = FALSE)'
RUN Rscript -e 'BiocManager::install("BSgenome", ask = FALSE)'

