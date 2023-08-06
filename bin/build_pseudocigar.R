library(argparse)

parser <- ArgumentParser(description='Create the pseudoCIGAR string using masked or unmasked ASVs')
parser$add_argument('--alignments', type="character", required=TRUE, help = "File containing aligned ASVs")

args <- parser$parse_args()
print(args)

library(stringr)
library(dplyr)
library(magrittr)

build_pseudoCIGAR_string <- function(reference, query) {

  # Transform strings into vectors of characters
  ref_chars <- str_split(reference, pattern = "")[[1]]
  query_chars <- str_split(query, pattern = "")[[1]]
  
  # Prepare a dataframe
  df <- tibble(position = seq_along(ref_chars),
               ref_char = ref_chars,
               query_char = query_chars)
  
  # Get positions that are in N-, NNNN, or N-N
  n_string <- str_locate_all(reference, pattern = "N+\\-+N*|N+")[[1]]
  n_vector <- unlist(mapply(seq, from = n_string[, "start"], to = n_string[, "end"]))

  # Indicate these positions are included in mask when building cigar
  df %<>% 
    mutate(masked = position %in% n_vector) %>%
    mutate(result = case_when(
      masked ~ paste0(n_string[, "start"], "N+", n_string[, "end"] - n_string[, "start"] + 1),
      ref_char != query_char & ref_char != "-" & query_char != "-" ~ paste0(position, query_char), 
      ref_char != query_char & ref_char == "-" & query_char != "-" & !position %in% ignore ~ paste0(position, "I=", query_char),
      ref_char != query_char & ref_char != "-" & query_char == "-" ~ paste0(position, "D=", ref_char),
      TRUE ~ ""
    ))
  
  # Combine the results into a single string
  cigar_str <- str_c(df$result, collapse = "")
  if (str_length(cigar_str) == 0) {
    cigar_str <- "." # if the cigar string is empty, replace with a dot
  }

  return (cigar_str)
}

df.aln = read.csv(args$alignments, sep="\t", header=T)
pseudo_cigar = df.aln %>%
  dplyr::mutate(pseudo_cigar = mapply(build_pseudoCIGAR_string, refseq, hapseq))

print(pseudo_cigar[1,])
write.table(pseudo_cigar,file="alignments.pseudocigar.txt",quote=F,sep="\t",col.names=T,row.names=F)
