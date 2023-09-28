library(logger)
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

load_library("dada2")
load_library("tidyverse")
load_library("argparse")

parser <- ArgumentParser(description='Create sequence table from DADA2 mergers')
parser$add_argument('--seqtab-paths', nargs='+', type="character", required=TRUE, help="Paths to RDS files containing DADA2 sequence tables")
parser$add_argument('--bimera-removal-method', type="character", required=TRUE, help="Method to remove bimeras")
parser$add_argument('--amplicon-file', type="character", required=TRUE, help="Path to the amplicon file")
parser$add_argument('--ncores', type="numeric", default=1, help="Number of threads to use")
parser$add_argument('--verbose', action='store_true', help="Add verbosity")
parser$add_argument('--dout', type="character", required=FALSE, help="Output directory")
parser$add_argument('--log-level', type="character", default = "INFO", help = "Log level. Default is INFO.")

args <- parser$parse_args()
# Set up logging
log_threshold(args$log_level)
log_appender(appender_console)
args_string <- paste(sapply(names(args), function(name) {
  paste(name, ":", args[[name]])
}), collapse = ", ")

log_debug(paste("Arguments parsed successfully:", args_string))

# Make a combined sequence table from sample specific sequence tables (RDS files)
seqtab_combined <- mergeSequenceTables(tables=args$seqtab_paths)

# Remove Bimeras - TODO parameterzie method
seqtab.nochim <- removeBimeraDenovo(seqtab_combined, method=args$bimera_removal_method, multithread=args$ncores, verbose=args$verbose)

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
