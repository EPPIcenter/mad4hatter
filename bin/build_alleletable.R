library(argparse)

parser <- ArgumentParser(description = "Create the allele table")
parser$add_argument("--amplicon-info", type = "character", required = TRUE, help = "File containing amplicon information")
parser$add_argument("--denoised-asvs", type = "character", required = TRUE, help = "File containing denoised ASVs")
parser$add_argument("--processed-asvs", type = "character", required = TRUE, help = "File containing processed ASVs")

args <- parser$parse_args()
print(args)

library(dplyr)

df.amplicon_info <- read.csv(args$amplicon_info, sep = "\t", header = T)
df.denoised <- read.csv(args$denoised_asvs, sep = "\t", header = T)
df.processed <- read.csv(args$processed_asvs, sep = "\t", header = T)

allele.table <- df.denoised %>%
  dplyr::inner_join(df.processed, by = c("sampleID", "locus" = "refid", "asv")) %>%
  select(sampleID, locus, asv, reads, pseudo_cigar) %>%
  group_by(sampleID, locus, asv) %>%
  mutate(reads = sum(reads)) %>%
  arrange(locus) %>%
  dplyr::rename(
    SampleID = sampleID,
    Locus = locus,
    ASV = asv,
    Reads = reads,
    PseudoCIGAR = pseudo_cigar
  ) %>%
  distinct()

# Filter amplicon info to only include 'amplicon' and 'pool' columns
df.amplicon_info_filtered <- df.amplicon_info %>%
  select(target_id, pool)
# Merge the amplicon info based on the locus and add the pool column
allele.table <- allele.table %>%
  dplyr::left_join(df.amplicon_info_filtered, by = c("Locus" = "target_id")) %>%
  dplyr::rename(Pool = pool) %>%
  select(SampleID, Locus, ASV, Reads, PseudoCIGAR, Pool)

write.table(allele.table, file = "allele_data.txt", quote = F, sep = "\t", col.names = T, row.names = F)
