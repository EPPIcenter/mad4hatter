library(Biostrings)
library(tidyverse)


files <- list.files(
  path = "Documents/repos/mad4hatter/",
  pattern = "^target_id_conversion_table\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)%>%
  # exclude any path containing /aad/
  .[!grepl("/aad/", .)]


combined <- map_dfr(files, function(f) {
  
  # Extract subfolder name (immediate parent folder)
  pool_name <- basename(dirname(f))
  
  read_tsv(f) %>%
    mutate(pool = pool_name)
}) %>% 
  group_by(new_name,old_name) %>% 
  summarize(pool = paste(pool,collapse = ",")) %>% 
  write_tsv("~/Documents/repos/mad4hatter/aad/target_id_conversion_table.tsv")



files <- list.files(
  path = "Documents/repos/mad4hatter/panel_information",
  pattern = "*_amplicon_info\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)%>%
  # exclude any path containing /aad/
  .[!grepl("/aad/", .)]



combined <- map_dfr(files, function(f) {
  
  # Extract subfolder name (immediate parent folder)
  pool_name <- basename(dirname(f))
  
  read_tsv(f) %>%
    mutate(pool = pool_name)
}) %>% 
  group_by(chrom,insert_start,insert_end,target_name,fwd_primer,rev_primer) %>% 
  summarize(pool = paste(pool,collapse = ",")) %>% 
  write_tsv("~/Documents/repos/mad4hatter/aad/amplicon_info.tsv")


files <- list.files(
  path = "Documents/repos/mad4hatter/panel_information",
  pattern = "*_reference\\.fasta$",
  recursive = TRUE,
  full.names = TRUE
)



# Read all FASTA sequences
all_seqs <- lapply(files, readDNAStringSet)  # or readAAStringSet for protein

# Combine all sequences
combined_seqs <- do.call(c, all_seqs)


# Keep only unique sequences and names
# unique() works on names + sequences together
# We'll use a data.frame approach to be safe
df <- data.frame(
  name = names(combined_seqs),
  seq  = as.character(combined_seqs),
  stringsAsFactors = FALSE
)

df_unique <- df[!duplicated(df$name) & !duplicated(df$seq), ]

# Convert back to DNAStringSet
unique_seqs <- DNAStringSet(df_unique$seq)
names(unique_seqs) <- df_unique$name

writeXStringSet(unique_seqs, "Documents/repos/mad4hatter/aad/references.fasta")

