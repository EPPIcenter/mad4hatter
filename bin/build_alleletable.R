library(argparse)

parser <- ArgumentParser(description = "Create the allele table")
parser$add_argument("--amplicon-info", type = "character", required = TRUE, help = "File containing amplicon information")
parser$add_argument("--denoised-asvs", type = "character", required = TRUE, help = "File containing denoised ASVs")
parser$add_argument("--processed-asvs", type = "character", required = TRUE, help = "File containing processed ASVs")
parser$add_argument("--processed-asvs-unmasked", type = "character", required = TRUE, help = "File containing processed ASVs without masking")
parser$add_argument("--aligned-asvs", type = "character", required = TRUE, help = "File containing aligned and masked ASVs")

args <- parser$parse_args()
print(args)

library(dplyr)

df.amplicon_info <- read.csv(args$amplicon_info, sep = "\t", header = T)
df.denoised <- read.csv(args$denoised_asvs, sep = "\t", header = T)
df.processed <- read.csv(args$processed_asvs, sep = "\t", header = T)
df.processed_unmasked <- read.csv(args$processed_asvs_unmasked, sep = "\t", header = T)
df.aligned_asvs <- read.csv(args$aligned_asvs, sep = "\t", header = T)

allele.table <- df.denoised %>%
  dplyr::inner_join(df.processed_unmasked, by = c("sampleID", "locus" = "refid", "asv")) %>%
  select(sampleID, locus, asv, reads, pseudo_cigar) %>%
  group_by(sampleID, locus, asv) %>%
  mutate(reads = sum(reads)) %>%
  arrange(locus) %>%
  dplyr::rename(
    experiment_sample_name = sampleID,
    target_id = locus,
    asv_raw = asv,
    reads = reads,
    pseudocigar_unmasked = pseudo_cigar
  ) %>%
  distinct()

allele.table <- allele.table %>%
  dplyr::inner_join(df.processed, by = c("experiment_sample_name"="sampleID", "target_id" = "refid", "asv_raw" = "asv")) %>%
  select(experiment_sample_name, target_id, asv_raw, reads, pseudocigar_unmasked, pseudo_cigar) %>%
  group_by(experiment_sample_name, target_id, asv_raw) %>%
  mutate(reads = sum(reads)) %>%
  arrange(target_id) %>%
  dplyr::rename(
    pseudocigar_masked = pseudo_cigar
  ) %>%
  distinct()

allele.table <- allele.table %>% 
  dplyr::inner_join(df.aligned_asvs, by = c("experiment_sample_name"="sampleID", "asv_raw" = "asv")) %>%
  select(experiment_sample_name, target_id, asv_raw, reads, pseudocigar_unmasked, pseudocigar_masked, masked_hapseq) %>%
  arrange(target_id) %>%
  dplyr::rename(
    asv_masked = masked_hapseq
  ) %>%
  distinct()

# Filter amplicon info to only include 'amplicon' and 'pool' columns
df.amplicon_info_filtered <- df.amplicon_info %>%
  select(target_id, pool)
# Merge the amplicon info based on the locus and add the pool column
allele.table <- allele.table %>%
  dplyr::left_join(df.amplicon_info_filtered, by = c("target_id" = "target_id")) %>%
  select(experiment_sample_name, target_id, asv_raw, pseudocigar_unmasked, asv_masked, pseudocigar_masked, pool, reads)

write.table(allele.table, file = "microhaplotypes_info.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

# Generate table in the old format
allele.table <- allele.table %>%
  select(experiment_sample_name, target_id, asv_raw, pseudocigar_masked, reads, pool) %>%
  dplyr::rename(
    SampleID = experiment_sample_name,
    Locus = target_id,
    ASV = asv_raw,
    Reads = reads,
    PseudoCIGAR = pseudocigar_masked,
    Pool = pool
)
write.table(allele.table, file = "allele_data.txt", quote = F, sep = "\t", col.names = T, row.names = F)
