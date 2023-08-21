library(argparse)
library(dplyr)
library(tidyr)
library(magrittr)

parser <- ArgumentParser(description='Calculate ASV Coverage')
parser$add_argument('--alleledata', type="character", required=TRUE, help="The allele table generated.")
parser$add_argument('--clusters', type="character", required=TRUE, help="DADA2 clusters")
parser$add_argument('--sample-coverage', type="character", help = "Sample coverage file from QC to append sample coverage statistics that are P. falciparum specific.")
parser$add_argument('--amplicon-coverage', type="character", help = "Amplicon coverage file from QC to append amplicon coverage statistics that are P. falciparum specific.")
parser$add_argument('--sample-coverage-out', type="character", help = "Postprocessed Sample coverage file", default = "sample_coverage_postprocessed.txt")
parser$add_argument('--amplicon-coverage-out', type="character", help = "Postprocessed Amplicon coverage file", default = "amplicon_coverage_postprocessed.txt")

args <- parser$parse_args()
print(args)

# setwd("~/Documents/GitHub/mad4hatter/work/e5/5b7adbb493c245f310327e0c87a340")

# args=list()
# args$alleledata="allele_data.txt"
# args$clusters="clusters.concatenated.collapsed.txt"
# args$sample_coverage="sample_coverage.txt"
# args$amplicon_coverage="amplicon_coverage.txt"
# args$sample_coverage_out="sample_coverage_postprocessed.txt"
# args$amplicon_coverage_out="amplicon_coverage_postprocessed.txt"


allele.data <- read.table(args$alleledata, header = TRUE, sep = "\t")
clusters <- read.table(args$clusters, header = TRUE, sep = "\t")


if (!is.null(args$sample_coverage) && file.exists(args$sample_coverage)) {

  sample.coverage = read.table(args$sample_coverage, header = TRUE, sep = "\t")
  print(head(sample.coverage))
  sample.coverage <- sample.coverage %>%
    pivot_wider(names_from = "X", values_from = "Reads")

  qc.postproc <- sample.coverage %>%
    left_join(clusters  %>%
      ungroup()  %>%
      select(sampleID,reads) %>%
      group_by(sampleID) %>%
      summarise(OutputDada2 = sum(reads)), by = c("SampleID" = "sampleID")
    )  %>%
    left_join(allele.data  %>%
      ungroup()  %>%
      select(sampleID,reads) %>%
      group_by(sampleID) %>%
      summarise(OutputPostprocessing = sum(reads)), by = c("SampleID" = "sampleID")
    )  %>%
    pivot_longer(cols = c(Input, `No Dimers`, Amplicons, OutputDada2, OutputPostprocessing))

  qc.postproc %<>% dplyr::rename("Stage" = "name", "Reads" = "value")


  write.table(qc.postproc, quote=F,sep='\t',col.names = TRUE, row.names = F, file = args$sample_coverage_out)
}

if (!is.null(args$amplicon_coverage) && file.exists(args$amplicon_coverage)) {
  amplicon.coverage <- read.table(args$amplicon_coverage, header = TRUE, sep = "\t")
  print(head(amplicon.coverage))

  qc.postproc <- amplicon.coverage  %>%
    left_join(clusters %>%
      group_by(sampleID,locus) %>%
      summarise(OutputDada2 = sum(reads)),
      by = c("SampleID" = "sampleID", "Locus" = "locus"),
      ) %>%
    left_join(allele.data %>%
      group_by(sampleID,locus) %>%
      summarise(OutputPostprocessing = sum(reads)),
          by = c("SampleID" = "sampleID", "Locus" = "locus"))
   qc.postproc$OutputDada2[is.na(qc.postproc$OutputDada2)] <- 0
   qc.postproc$OutputPostprocessing[is.na(qc.postproc$OutputPostprocessing)] <- 0

  write.table(qc.postproc, quote=F,sep='\t',col.names = TRUE, row.names = F, file = args$amplicon_coverage_out)
}
