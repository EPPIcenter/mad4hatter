library(argparse)

parser <- ArgumentParser(description='Post process dada2 inferred sequences')
parser$add_argument('--alignments', type="character", required=TRUE, help="RDS Clusters from DADA2. This is the main output from the DADA module.")
parser$add_argument('--masks', type="character", nargs="+")
parser$add_argument('--parallel', action='store_true')
parser$add_argument('--n-cores', type = 'integer', default = -1, help = "Number of cores to use. Ignored if running parallel flag is unset.")

args <- parser$parse_args()
print(args)

library(Biostrings)
library(dplyr)
library(purrr)
library(doMC)
library(parallel)
library(stringr)
library(magrittr)

# setwd("/home/bpalmer/Documents/GitHub/mad4hatter/work/ac/bc2fb8f020f492c7de5baa45f44fb7")
# args=list()
# masks=""

df_aln=read.table(args$alignments, header=T)

## Setup parallel backend if asked

if (args$parallel) {
  n_cores <- ifelse(args$n_cores <= 0, detectCores(), args$n_cores)
  registerDoMC(n_cores)
} else {
  registerDoSEQ()
}

# Function to find 'N' positions
get_masked_pos <- function(sequence) {
  out = unlist(gregexpr("N", sequence, ignore.case = TRUE)) # Using 'gregexpr' to find the positions
  if (length(out) == 1 && out == -1) {
    return (NA)
  } else {
    return (out)
  }
}

fasta_positions <- lapply(args$masks, function(fasta_path) {
  fasta <- readDNAStringSet(fasta_path)
  # Create a data frame for each sequence
  positions_df <- data.frame(refid=names(fasta), stringsAsFactors=FALSE)
  positions_df$pos <- lapply(fasta, function(x) get_masked_pos(as.character(x))) # Store the positions as a list within the dataframe
  names(positions_df)[2] <- fasta_path # Name the second column as the fasta file path
  positions_df
})

# Merge all the data frames
merged <- Reduce(function(df1, df2) merge(df1, df2, by="refid", all=TRUE), fasta_positions)
merged$all_positions <- apply(merged[-which(names(merged) %in% "refid")], 1, function(x) unlist(x[!is.na(unlist(x))]))
merged = merged[, c("refid", "all_positions")]

# Now let's mask the aligned reference sequences
# first make a data frame with all the unique aligned sequences and join with the data frame that has the masked references
df_aln_references = df_aln  %>%
  select(refid,refseq) %>%
  distinct() %>%
  left_join(merged,by="refid")


mask_aligned_sequence = function(pos, refseq) {

  # return early if there is no need to mask
  if (is.null(pos) || (is.list(pos) && length(pos) == 0)) {
    return (refseq)
  }

  pattern <- gregexpr("[^-]+|-+", refseq)
  split_string <- regmatches(refseq, pattern)[[1]]
  string_homo_tr_ref_seq = str_remove_all(refseq, "-")

  # build the reference sequence using collapsed positions
  string_homo_tr_ref_seq = strsplit(as.character(string_homo_tr_ref_seq), "")[[1]]
  string_homo_tr_ref_seq[unlist(pos)] <- "N"
  string_homo_tr_ref_seq <- paste(string_homo_tr_ref_seq, collapse = "")

  newstring = character()

  for(chunk in split_string){
    if(grepl("-",chunk)){
      newstring=paste0(newstring,chunk)
    }else{
      newstring=paste0(newstring,substr(string_homo_tr_ref_seq,1,nchar(chunk)))
      string_homo_tr_ref_seq = substr(string_homo_tr_ref_seq,nchar(chunk)+1,nchar(string_homo_tr_ref_seq))
    }
  }
  return(newstring)
}

remove_dashes_between_N = function(newstring){
  NdashN = unique(regmatches(newstring,gregexpr("N-+N", newstring))[[1]])
  newstring_nodash = newstring
  for(i in NdashN){
    newstring_nodash=gsub(i,strrep("N",nchar(i)),newstring_nodash)
  }
  return(newstring_nodash)
  }

pos=df_aln_references[3,]$all_positions
refseq=df_aln_references[3,]$refseq

df_aln_references %<>%
  mutate(masked_aligned_refseq = mapply(mask_aligned_sequence,all_positions,refseq)) %>%
  mutate(masked_aligned_nodashinN_refseq = sapply(masked_aligned_refseq,remove_dashes_between_N))

#update df_aln
df_aln %<>% left_join(df_aln_references,by=c("refid","refseq")) %>%
  select(sampleID, asv, refid, hapseq, masked_aligned_nodashinN_refseq, score, indels) %>%
  dplyr::rename(refseq = masked_aligned_nodashinN_refseq)

write.table(df_aln,file="masked.alignments.txt",quote=F,sep="\t",col.names=T,row.names=F)
