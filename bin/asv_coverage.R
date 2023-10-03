library(logger)
log_threshold(WARN)
log_appender(appender_console)

load_library <- function(library_name) {
  output <- capture.output({
    suppressWarnings({
      library(library_name, character.only = TRUE)
    })
  }, type = "message")
  
  # Separate warnings from messages
  warnings <- warnings()
  
  # Log messages
  if(length(output) > 0) {
    log_info(paste("Message from", library_name, ":", paste(output, collapse = "; ")))
  }
  
  # Log warnings
  if(length(warnings) > 0) {
    log_warn(paste("Warning from", library_name, ":", paste(warnings, collapse = "; ")))
  }
}

load_library("argparse")
load_library("dplyr")
load_library("tidyr")
load_library("magrittr")


parser <- ArgumentParser(description='Calculate ASV Coverage')
parser$add_argument('--alleledata', type="character", required=TRUE, help="The allele table generated.")
parser$add_argument('--clusters', type="character", required=TRUE, help="DADA2 clusters")
parser$add_argument('--sample-coverage', type="character", help = "Sample coverage file from QC to append sample coverage statistics that are P. falciparum specific.")
parser$add_argument('--amplicon-coverage', type="character", help = "Amplicon coverage file from QC to append amplicon coverage statistics that are P. falciparum specific.")
parser$add_argument('--sample-coverage-out', type="character", help = "Postprocessed Sample coverage file", default = "sample_coverage_postprocessed.txt")
parser$add_argument('--amplicon-coverage-out', type="character", help = "Postprocessed Amplicon coverage file", default = "amplicon_coverage_postprocessed.txt")
parser$add_argument('--log-level', type="character", default = "INFO", help = "Log level. Default is INFO.")

args <- parser$parse_args()

log_level_arg <- match.arg(args$log_level, c("DEBUG", "INFO", "WARN", "ERROR", "FATAL"))
log_threshold(log_level_arg)

args_string <- paste(sapply(names(args), function(name) {
  paste(name, ":", args[[name]])
}), collapse = ", ")

log_debug(paste("Arguments parsed successfully:", args_string))
log_info("Reading allele data from {args$alleledata}")
allele.data <- read.table(args$alleledata, header = TRUE, sep = "\t")

log_info("Reading clusters from {args$clusters}")
clusters <- read.table(args$clusters, header = TRUE, sep = "\t")


if (!is.null(args$sample_coverage) && file.exists(args$sample_coverage)) {
  log_info("Processing sample coverage data...")
  sample.coverage = read.table(args$sample_coverage, header = TRUE, sep = "\t")
  print(str(sample.coverage))
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
    mutate(across(everything(), as.character)) %>%    
    pivot_longer(cols = c(Input, `No Dimers`, Amplicons, OutputDada2, OutputPostprocessing))

  qc.postproc %<>% dplyr::rename("Stage" = "name", "Reads" = "value")


  write.table(qc.postproc, quote=F,sep='\t',col.names = TRUE, row.names = F, file = args$sample_coverage_out)
  log_info("Writing postprocessed data to {args$sample_coverage_out}")
} else {
  log_error("Sample coverage file not found or not specified.")
}

if (!is.null(args$amplicon_coverage) && file.exists(args$amplicon_coverage)) {
  amplicon.coverage <- read.table(args$amplicon_coverage, header = TRUE, sep = "\t")
  log_info("Processing amplicon coverage data...")
  log_debug(str(amplicon.coverage))

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
  log_info("Writing postprocessed data to {args$amplicon_coverage_out}")
} else {
  log_error("Amplicon coverage file not found or not specified.")
}
