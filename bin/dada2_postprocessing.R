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

load_library("lobstr")
# Define profiling function
profile_function <- function(func, ...) {
  # Record the memory size of the arguments
  mem_args <- lobstr::obj_size(...)

  # Variable to store messages
  messages_captured <- character(0)

  # Record start time and memory
  time_start <- proc.time()[["elapsed"]]
  mem_start <- lobstr::mem_used()

  # Evaluate the function while capturing its output and messages
  output_captured <- capture.output({
      withCallingHandlers({
          result <- func(...)
      }, message = function(m) {
          # Capture just the content of the message
          messages_captured <<- c(messages_captured, conditionMessage(m))
      })
  })

  # Record end time and memory
  time_end <- proc.time()[["elapsed"]]
  mem_end <- lobstr::mem_used()

  # Combine standard output and messages
  all_captured <- c(output_captured, messages_captured)

  # Return results including captured output and messages
  list(
      result = result, 
      output = all_captured,
      mem_args = mem_args,
      mem_diff = (mem_end - mem_start) - mem_args, 
      duration = time_end - time_start
  )
}

load_library("dada2")
load_library("tidyverse")
load_library("argparse")

parser <- ArgumentParser(description='Create sequence table from DADA2 mergers')
parser$add_argument('--seqtab-paths', nargs='+', type="character", required=TRUE, help="Paths to RDS files containing DADA2 sequence tables")
parser$add_argument('--bimera-removal-method', type="character", required=TRUE, help="Method to remove bimeras")
parser$add_argument('--amplicon-file', type="character", required=TRUE, help="Path to the amplicon file")
parser$add_argument('--ncores', type="numeric", default=1, help="Number of threads to use")
parser$add_argument('--dout', type="character", required=FALSE, help="Output directory")
parser$add_argument('--log-level', type="character", default = "INFO", help = "Log level. Default is INFO.")

args <- parser$parse_args()
# Set up logging
log_level_arg <- match.arg(args$log_level, c("DEBUG", "INFO", "WARN", "ERROR", "FATAL"))
log_threshold(log_level_arg)

args_string <- paste(sapply(names(args), function(name) {
  paste(name, ":", args[[name]])
}), collapse = ", ")

log_debug(paste("Arguments parsed successfully:", args_string))

# Profile and capture output for mergeSequenceTables
merge_profile <- profile_function(function(tables) {
    mergeSequenceTables(tables=tables)
}, args$seqtab_paths)

# Log profiling results
log_info("mergeSequenceTables duration: {merge_profile$duration}")
log_info("Memory change during mergeSequenceTables: {merge_profile$mem_diff}")
log_info("Memory used by arguments to mergeSequenceTables: {merge_profile$mem_args}")

# Use the result from the profiled function
seqtab_combined <- merge_profile$result
log_info("{length(args$seqtab_paths)} sequence tables merged.")

# Profile and capture output for removeBimeraDenovo
bimera_profile <- profile_function(function(seqtab_combined) {
    removeBimeraDenovo(seqtab_combined, method=args$bimera_removal_method, multithread=args$ncores, verbose=TRUE)
}, seqtab_combined)

# Log profiling results
log_info("removeBimeraDenovo duration: {bimera_profile$duration} seconds")
log_info("Memory change during removeBimeraDenovo: {bimera_profile$mem_diff}")

# Use the result from the profiled function
seqtab.nochim <- bimera_profile$result

# Log the captured output
if(length(bimera_profile$output) > 0) {
    log_info("Output from removeBimeraDenovo:\n{paste(bimera_profile$output, collapse = '\n')}")
}

### Create clusters output file with quality control
seqtab.nochim.df = as.data.frame(seqtab.nochim)
seqtab.nochim.df$sample = rownames(seqtab.nochim)
seqtab.nochim.df[seqtab.nochim.df==0]=NA

amplicon.table=read.table(args$amplicon_file, header = TRUE, sep = "\t")

# find amplicons (use to select)
pat=paste(amplicon.table$amplicon, collapse="|") # find amplicons specified in amplicon table

# create regex to extract sampleID (to be used with remove)
pat.sampleID=paste(sprintf("^%s_", amplicon.table$amplicon), collapse="|") # find amplicons specified in amplicon table (with _)
pat.sampleID=paste(c(pat.sampleID, "_trimmed_merged.fastq.gz$"),collapse="|")  # remove _trimmed from end

seqtab.nochim.df = seqtab.nochim.df %>%
  pivot_longer(cols = seq(1,ncol(seqtab.nochim)),names_to = "asv",values_to = "reads",values_drop_na=TRUE) %>%
  mutate(locus = str_extract(sample, pat)) %>%
  mutate(sampleID = str_remove_all(sample, pat.sampleID)) %>%
  select(sampleID,locus,asv,reads)

temp = seqtab.nochim.df %>% select(locus,asv) %>% distinct()
loci =unique(temp$locus)
k=1
allele.sequences = data.frame(locus = seq(1,nrow(temp)),allele = seq(1,nrow(temp)),sequence = seq(1,nrow(temp)))
for(i in seq(1,length(loci))){
  temp2 = temp %>% filter(locus==loci[i])
  for(j in seq(1,nrow(temp2))){
    allele.sequences$locus[k+j-1] = loci[i]
    allele.sequences$allele[k+j-1] = paste0(loci[i],".",j)
    allele.sequences$sequence[k+j-1] = temp2$asv[j]
  }
  k=k+nrow(temp2)
}

allele.data = seqtab.nochim.df %>%
  mutate(sampleID = str_remove_all(sampleID, pat = "_trimmed")) %>%
  left_join(allele.sequences %>% select(-locus),by=c("asv"="sequence")) %>%
  group_by(sampleID,locus) %>%
  mutate(norm.reads.locus = reads/sum(reads))%>%
  mutate(n.alleles = n()) %>%
  ungroup()

output_dir <- if(is.null(args$dout)) {
  "."
} else {
  args$dout
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive=TRUE)
}

write.table(allele.data,file=file.path(output_dir, "dada2.clusters.txt") ,quote=F,sep="\t",col.names=T,row.names=F)
