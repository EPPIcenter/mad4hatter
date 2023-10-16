library(argparse)

parser <- ArgumentParser(description='Create the allele table')
parser$add_argument('--denoised-asvs', type="character", required=TRUE, help = "File containing denoised ASVs")
parser$add_argument('--processed-asvs', type="character", required=TRUE, help = "File containing processed ASVs")

args <- parser$parse_args()
print(args)

library(dplyr)
# setwd("/home/bpalmer/Documents/GitHub/mad4hatter/work/50/8c3ac24d1982d59265e1efcf64488e")
# args=list()
# args$denoised_asvs = "clusters.concatenated.collapsed.txt"
# args$processed_asvs = "alignments.pseudocigar.txt"

df.denoised = read.csv(args$denoised_asvs, sep="\t", header=T)
df.processed = read.csv(args$processed_asvs, sep="\t", header=T)

allele.table = df.denoised %>%
  dplyr::inner_join(df.processed, by=c("sampleID", "locus"="refid", "asv")) %>%
  select(sampleID, locus, asv, reads, allele, pseudo_cigar) %>%
  group_by(sampleID, locus, asv) %>%
  mutate(reads=sum(reads)) %>%
  arrange(locus) %>%
  dplyr::rename(
    SampleID = sampleID,
    Locus = locus,
    ASV = asv,
    Reads = reads,
    Allele = allele,
    PseudoCIGAR = pseudo_cigar
  ) %>%
  distinct()

write.table(allele.table,file="allele_data.txt",quote=F,sep="\t",col.names=T,row.names=F)
