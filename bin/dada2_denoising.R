library(dada2)
library(tidyverse)
library(argparse)

parser <- ArgumentParser(description='DADA2 Denoising')
parser$add_argument('--derep-1', type="character", required=TRUE, nargs="+", help="Path to RDS file containing dereplicated forward sequences")
parser$add_argument('--derep-2', type="character", required=TRUE, nargs="+", help="Path to RDS file containing dereplicated reverse sequences")
parser$add_argument('--error-model-1', type="character", required=TRUE, help="Path to the forward error model RDS file")
parser$add_argument('--error-model-2', type="character", required=TRUE, help="Path to the reverse error model RDS file")
parser$add_argument('--ncores', type="numeric", required=TRUE, help="Number of threads to use")
parser$add_argument('--pool', type="character", default="false")
parser$add_argument('--band-size', type='integer', default=16)
parser$add_argument('--omega-a', type='double', default=1e-120)
parser$add_argument('--just-concatenate', action='store_true')
parser$add_argument('--use-quals', type="character", default="false")
parser$add_argument('--maxEE', type="integer", default=2)
parser$add_argument('--self-consist', action='store_true')
parser$add_argument('--omega-c', type='double', default=1e-40)
parser$add_argument('--ampliconFILE', type="character", required=TRUE, help="Path to amplicon info file")
parser$add_argument('--verbose', action="store_true", help="Verbose")

args <- parser$parse_args()
print(args)

# Load error models
err_model_F <- readRDS(args$error_model_1)
err_model_R <- readRDS(args$error_model_2)

# Function to gather derep objects from a list of files
gather_dereps <- function(paths) {
    all_dereps <- list()
    for(path in paths) {
        derep <- readRDS(path)
        # Extract the sample name using a regex pattern
        sample_name <- sub("(.*)_trimmed.*", "\\1", basename(path))
        all_dereps[[sample_name]] <- derep
    }
    all_dereps
}

# Gather all derep objects
derepFs <- gather_dereps(args$derep_1)
derepRs <- gather_dereps(args$derep_2)

# DADA2 denoising
dadaFs <- dada(derepFs, err=err_model_F, multithread=args$ncores)
dadaRs <- dada(derepRs, err=err_model_R, multithread=args$ncores)

# Save the DADA2 denoised data
save(dadaFs, file="dada_denoised_F.RData")
save(dadaRs, file="dada_denoised_R.RData")

# Filter samples based on the provided condition
filter_samples <- function(samples, condition) {
  samples[names(samples) %in% rownames(condition)]
}

get_sample_name <- function(filename, loci_list) {
  # Remove locus from the beginning
  no_locus_filename <- gsub(paste0("^(", paste(loci_list, collapse="|"), ")_"), "", filename)
  
  # Remove everything past '_trimmed'
  sample_name <- strsplit(no_locus_filename, "_trimmed")[[1]][1]
  
  return(sample_name)
}

merge_with_concatenation <- function(dadaFs, derepFs, dadaRs, derepRs, amplicon_info) {
  # Perform filtering
  overlap_condition <- amplicon_info %>% filter((ampInsert_length + 10) < sum.mean.length.reads)
  no_overlap_condition <- amplicon_info %>% filter((ampInsert_length + 10) >= sum.mean.length.reads)
  
  # Merge for overlap and no overlap conditions
  mergers.overlap <- perform_merging(filter_samples(dadaFs, overlap_condition),
                                     filter_samples(derepFs, overlap_condition),
                                     filter_samples(dadaRs, overlap_condition),
                                     filter_samples(derepRs, overlap_condition),
                                     justConcatenate=FALSE)

  mergers.no.overlap <- perform_merging(filter_samples(dadaFs, no_overlap_condition),
                                        filter_samples(derepFs, no_overlap_condition),
                                        filter_samples(dadaRs, no_overlap_condition),
                                        filter_samples(derepRs, no_overlap_condition),
                                        justConcatenate=TRUE)

  # Combine the results
  c(mergers.overlap, mergers.no.overlap)[names(dadaFs)]
}

perform_merging <- function(dadaFs_selected, derepFs_selected, dadaRs_selected, derepRs_selected, justConcatenate) {
  mergePairs(dadaFs_selected, derepFs_selected, dadaRs_selected, derepRs_selected, 
             verbose = TRUE, 
             justConcatenate = justConcatenate, 
             trimOverhang = TRUE,
             minOverlap = 10,
             maxMismatch = 1)
}

amplicon_data <- read.table(args$ampliconFILE, header=T)
loci_list <- amplicon_data$amplicon

if (args$just_concatenate) {
  amplicon.info = data.frame(
    names = sapply(strsplit(names(dadaFs),"_trimmed"),"[",1)
  )

  # Create amplicon info
  amplicon.info <- amplicon.info %>%
    mutate(names=sapply(str_split(names,'_S(\\d+)'),head,1)) %>% 
    mutate(amplicon=unlist(lapply(str_split(names,'_'), function(x) { paste(x[1:3], collapse = "_") }))) %>%
    inner_join(amplicon_data %>%
               select(amplicon, ampInsert_length), by = c("amplicon")) %>%
    select(amplicon, ampInsert_length) %>%
    mutate(
      sum.mean.length.reads = sapply(sapply(unlist(lapply(dadaFs,"[","sequence"),recursive = F),nchar),mean)+
        sapply(sapply(unlist(lapply(dadaRs,"[","sequence"),recursive = F),nchar),mean)
    )

  rownames(amplicon.info) = names(dadaFs)

  mergers <- merge_with_concatenation(dadaFs, derepFs, dadaRs, derepRs, amplicon.info)
  
  sample_name <- get_sample_name(names(derepFs)[1], loci_list)  # extract sample name
  saveRDS(mergers, file = paste0(sample_name, "_merged.RDS"))
} else {
  # Perform normal merging
  mergers <- perform_merging(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=FALSE)
  
  sample_name <- get_sample_name(names(derepFs)[1], loci_list)  # extract sample name
  saveRDS(mergers, file = paste0(sample_name, "_merged.RDS"))
}

seqtab <- makeSequenceTable(mergers)

sample_name <- get_sample_name(names(derepFs)[1], loci_list)  # extract sample name
saveRDS(seqtab, file = paste0(sample_name, "_seqtab.RDS"))
