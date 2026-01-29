library(argparse)
library(magrittr)

parser <- ArgumentParser(description = "Create the allele table")
parser$add_argument("--amplicon-info", type = "character", required = TRUE, help = "File containing amplicon information")
parser$add_argument("--denoised-asvs", type = "character", required = TRUE, help = "File containing denoised ASVs")
parser$add_argument("--masked-pseudocigar-table", type = "character", required = TRUE, help = "File containing processed ASVs")
parser$add_argument("--unmasked-pseudocigar-table", type = "character", required = TRUE, help = "File containing unmasked ASVs")
parser$add_argument("--masked-asv-table", type = "character", required = TRUE, help = "File containing aligned and masked ASVs")

args <- parser$parse_args()
print(args)

library(dplyr)

df.amplicon_info <- read.csv(args$amplicon_info, sep = "\t", header = T)
df.denoised <- read.csv(args$denoised_asvs, sep = "\t", header = T)
df.masked_pseudocigar <- read.csv(args$masked_pseudocigar_table, sep = "\t", header = T)
df.unmasked_pseudocigar <- read.csv(args$unmasked_pseudocigar_table, sep = "\t", header = T)
df.masked_asv_table <- read.csv(args$masked_asv_table, sep = "\t", header = T)

df.masked_pseudocigar %<>%
  dplyr::rename(pseudocigar_masked = pseudo_cigar)

# TODO: Add in a catch that if masked_hapseq is not in the aligned_asvs table then use the unmasked hapseq
allele.table <- df.denoised %>%
  dplyr::inner_join(df.unmasked_pseudocigar, by = c("sampleID", "locus" = "refid", "asv")) %>%
  dplyr::inner_join(df.masked_pseudocigar, by = c("sampleID", "locus" = "refid", "asv")) %>%
  dplyr::inner_join(df.masked_asv_table, by = c("sampleID", "locus" = "refid", "asv")) %>%
  select(sampleID, locus, asv, reads, pseudo_cigar, pseudocigar_masked, masked_hapseq) %>%
  #group_by(sampleID, locus, asv) %>%
  #mutate(reads = sum(reads)) %>%
  arrange(locus) %>%
  dplyr::rename(
    sample_name = sampleID,
    target_name = locus,
    pseudocigar_unmasked = pseudo_cigar,
    asv_masked = masked_hapseq
  ) %>%
  distinct()

# Filter amplicon info to only include 'amplicon' and 'pool' columns
df.amplicon_info_filtered <- df.amplicon_info %>%
  select(target_id, pool)
# Merge the amplicon info based on the locus and add the pool column
allele.table <- allele.table %>%
  dplyr::left_join(df.amplicon_info_filtered, by = c("target_name" = "target_id")) %>%
  select(sample_name, target_name, asv, pseudocigar_unmasked, asv_masked, pseudocigar_masked, reads, pool)
write.table(allele.table, file = "allele_data.txt", quote = F, sep = "\t", col.names = T, row.names = F)

allele.table.collapsed <- allele.table %>%
  select(sample_name, target_name, asv_masked, pseudocigar_masked, reads, pool) %>%
  group_by(sample_name, target_name, asv_masked, pseudocigar_masked) %>%
  mutate(reads = sum(reads)) %>%
  arrange(target_name) %>%
  distinct()

write.table(allele.table.collapsed, file = "allele_data_collapsed.txt", quote = F, sep = "\t", col.names = T, row.names = F)