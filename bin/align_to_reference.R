library(argparse)

parser <- ArgumentParser(description='Aligns sequences (DADA2 clusters) against a reference of known sequences')
parser$add_argument('--clusters', type="character", required=TRUE, help="RDS Clusters from DADA2. This is the main output from the DADA module.")
parser$add_argument('--refseq-fasta', type="character")
parser$add_argument('--alignment-threshold', type="integer", default = 60)
parser$add_argument('--n-cores', type = 'integer', default = 1, help = "Number of cores to use. Ignored if running parallel flag is unset.")
parser$add_argument('--amplicon-table', type="character", required=TRUE, help = "Amplicon table with primer pools. This is used to organize the sequence table by amplicon.")

args <- parser$parse_args()
print(args)

library(stringr)
library(dplyr)
library(foreach)
library(parallel)
library(BSgenome)
library(tidyr)
library(doMC)
library(tibble)
library(ggplot2)
library(Biostrings)
library(magrittr)

# register number of cores to use
registerDoMC(args$n_cores)

# read the sequences from the reference genome (already extracted into a fasta file)
clusters <- read.table(args$clusters, header = TRUE)
ref_sequences <- readDNAStringSet(args$refseq_fasta)

# DADA2 should have corrected any substitution errors during sequencing. Therefore, any 
# mismatches are expected to be real. A small penalty is added when there is a mismatching base
# to catch large differences that are likely not real. Gap penalties are also set below to 
# help filter sequences that are truly different from the reference. 
sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)

# Perform alignment by unique ASV by locus
unique_clusters <- clusters %>%
  distinct(locus, asv)

df_aln <- foreach(ii = 1:nrow(unique_clusters), .combine = "bind_rows") %dopar% {
  # the alignment is performed only with the reference sequence for the corresponding locus as we have this information from the demultiplexing step
  refseq <- ref_sequences[unique_clusters$locus[ii]]
  aln <- pairwiseAlignment(refseq, str_remove_all(unique_clusters$asv[ii], "N"), substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE)
  patt <- c(alignedPattern(aln), alignedSubject(aln))
  ind <- sum(str_count(as.character(patt), "-"))
  data.frame(
    locus = unique_clusters$locus[ii],
    asv = unique_clusters$asv[ii],
    hapseq = as.character(patt)[2],
    refseq = as.character(patt)[1],
    score = score(aln),
    indels = ind
  )
}

# Merge the alignment results back to the original clusters data
df_aln_merged <- clusters %>%
  dplyr::left_join(df_aln, by = c("locus", "asv")) %>%
  dplyr::rename(refid = locus) %>%
  dplyr::select(sampleID, asv, hapseq, refseq, refid, score, indels)

write.table(df_aln_merged, file = "alignments.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
