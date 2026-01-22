library(argparse)
library(magrittr)

parser <- ArgumentParser(description = "Create the allele table")
parser$add_argument("--amplicon-info", type = "character", required = TRUE, help = "File containing amplicon information")
parser$add_argument("--denoised-asvs", type = "character", required = TRUE, help = "File containing denoised ASVs")
parser$add_argument("--processed-asvs-masked", type = "character", required = TRUE, help = "File containing processed ASVs")
parser$add_argument("--processed-asvs-unmasked", type = "character", required = TRUE, help = "File containing unmasked ASVs")
parser$add_argument("--aligned-asvs", type = "character", required = TRUE, help = "File containing aligned and masked ASVs")

args <- parser$parse_args()
print(args)

library(dplyr)

df.amplicon_info <- read.csv(args$amplicon_info, sep = "\t", header = T)
df.denoised <- read.csv(args$denoised_asvs, sep = "\t", header = T)
df.processed_masked <- read.csv(args$processed_asvs_masked, sep = "\t", header = T)
df.processed_unmasked <- read.csv(args$processed_asvs_unmasked, sep = "\t", header = T)
df.aligned_asvs <- read.csv(args$aligned_asvs, sep = "\t", header = T)

df.processed_masked %<>%
  dplyr::rename(pseudocigar_masked = pseudo_cigar)

# TODO: Add in a catch that if masked_hapseq is not in the aligned_asvs table then use the unmasked hapseq
allele.table <- df.denoised %>%
  dplyr::inner_join(df.processed_unmasked, by = c("sampleID", "locus" = "refid", "asv")) %>%
  dplyr::inner_join(df.processed_masked, by = c("sampleID", "locus" = "refid", "asv")) %>%
  dplyr::inner_join(df.aligned_asvs, by = c("sampleID", "locus" = "refid", "asv")) %>%
  select(sampleID, locus, asv, reads, pseudo_cigar, pseudocigar_masked, masked_hapseq) %>%
  group_by(sampleID, locus, asv) %>%
  mutate(reads = sum(reads)) %>%
  arrange(locus) %>%
  dplyr::rename(
    sample_name = sampleID,
    target_name = locus,
    asv_raw = asv,
    read_count = reads,
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
  select(sample_name, target_name, asv_raw, pseudocigar_unmasked, asv_masked, pseudocigar_masked, read_count, pool)
write.table(allele.table, file = "microhaplotypes_info.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

# Convert to the old format
allele.table_old <- allele.table %>%
  select(sample_name, target_name, asv_raw, read_count, pseudocigar_masked, pool) %>%
  dplyr::rename(SampleID = sample_name, Locus = target_name, ASV = asv_raw, PseudoCIGAR = pseudocigar_masked, Reads = read_count, Pool = pool)

write.table(allele.table_old, file = "allele_data.txt", quote = F, sep = "\t", col.names = T, row.names = F)
