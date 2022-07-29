library(moire)
library(tidyverse)
args = commandArgs(trailingOnly=T)

haps <- readRDS(args[1])

haps2 = haps %>% select(sampleID,locus) %>% distinct()


haps.targets = haps2 %>% group_by(sampleID) %>%  dplyr::summarise(n = n()) 

#samples_keep = haps.targets %>% filter(n>90) %>% select(s_Sample)


haps.clean = haps %>% select(sampleID,locus,allele,reads) %>% mutate_at("reads",as.numeric)

data.moire.table = haps.clean %>% select(sampleID,locus,allele) %>% ungroup()
colnames(data.moire.table) = c("sample_id", "locus", "allele")


temp = data.moire.table  %>% dplyr::select(-sample_id) %>% distinct() %>% group_by(locus) %>% dplyr::summarise(n = n()) %>% filter(n==1)
data.moire.table = data.moire.table %>% filter(!(locus %in% temp$locus ))

data.moire = load_long_form_data(data.moire.table)

results.moire = run_mcmc(
  data.moire$data,data.moire$sample_ids,data.moire$loci,data.moire$is_missing,
  verbose = T, burnin = 1e4, samples = 1e4, thin = 1,
  eps_pos_alpha = 1, eps_pos_beta = 99,
  eps_neg_alpha = 10, eps_neg_beta = 990, allele_freq_vars = 1
)


# Estimate the COI for each sample
coi_summary <- moire::summarize_coi(results.moire)

# We can also summarize statistics about the allele frequency distribution
he_summary <- moire::summarize_he(results.moire)
allele_freq_summary <- moire::summarize_allele_freqs(results.moire)

  write.table(coi_summary,file="coi_summary.txt",quote=F,sep="\t",col.names=T,row.names=F)
