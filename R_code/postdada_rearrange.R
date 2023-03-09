library(argparse)

parser <- ArgumentParser(description='Post process dada2 inferred sequences')
parser$add_argument('--homopolymer-threshold', type="integer",
                   help='homopolymer threshold to begin masking')
parser$add_argument('--refseq-fasta', type="character")
parser$add_argument('--masked-fasta', type="character")
parser$add_argument('--dada2-output', type="character", required = TRUE)
parser$add_argument('--alignment-threshold', type="integer", default = 60)
parser$add_argument('--parallel', action='store_true')
parser$add_argument('--n-cores', type = 'integer', default = -1)

args <- parser$parse_args()
print(args)

library(stringr)
library(dplyr)
library(dada2)
library(foreach)
library(parallel)
library(muscle)
library(BSgenome)
library(tidyr)
library(doMC)
library(tibble)


# FOR DEBUGGING
# args <- list()
# args$homopolymer_threshold <- 5
# args$refseq_fasta <- "v4_refseq.fasta"
# args$masked_fasta <- "v4_refseq.fasta.2.7.7.80.10.25.3.mask"
# args$dada2_output <- "seqtab.nochim.RDS"
# args$parallel <- FALSE
# args$alignment_threshold <- 60

## Postprocessing QC

seqtab.nochim <- readRDS(args$dada2_output)
seqtab.nochim.df = as.data.frame(seqtab.nochim)
seqtab.nochim.df$sample = rownames(seqtab.nochim)
seqtab.nochim.df[seqtab.nochim.df==0]=NA
pat="-1A_|-1B_|-1_|-2_|-1AB_|-1B2_"
seqtab.nochim.df = seqtab.nochim.df %>%
  pivot_longer(cols = seq(1,ncol(seqtab.nochim)),names_to = "asv",values_to = "reads",values_drop_na=TRUE) %>%
  mutate(locus = paste0(sapply(strsplit(sample,"_"),"[",1),"_",sapply(strsplit(sample,"_"),"[",2),"_",sapply(strsplit(sample,"_"),"[",3)))%>%
  mutate(sampleID = sapply(strsplit(sapply(strsplit(sample,pat),"[",2),"_trimmed"),"[",1)) %>%
  select(sampleID,locus,asv,reads)

temp = seqtab.nochim.df %>% select(locus,asv) %>% distinct()
loci =unique(temp$locus)
k=1
allele.sequences = data.frame(locus = seq(1,nrow(temp)),allele = seq(1,nrow(temp)),sequence = seq(1,nrow(temp)))
for(i in seq(1,length(loci))){
  temp2 = temp %>% filter(locus==loci[i])
  for(j in seq(1,nrow(temp2))){
    allele.sequences$locus[k+j-1] = loci[i]
    allele.sequences$allele[k+j-1] = paste0(loci[i],".",j)
    allele.sequences$sequence[k+j-1] = temp2$asv[j]
  }
  k=k+nrow(temp2)
}

allele.data = seqtab.nochim.df %>%
  left_join(allele.sequences %>% select(-locus),by=c("asv"="sequence")) %>%
  group_by(sampleID,locus,allele) %>%
  mutate(norm.reads.allele = reads/sum(reads))%>%
  group_by(sampleID,locus) %>%
  mutate(norm.reads.locus = reads/sum(reads))%>%
  mutate(n.alleles = n())

saveRDS(allele.data,file="pre_processed_allele_table.RDS")
write.table(allele.data,file="pre_processed_allele_table.txt",quote=F,sep="\t",col.names=T,row.names=F)

seqtab.nochim <- readRDS(args$dada2_output)

## I. Check for non overlapping sequences.

sequences = colnames(seqtab.nochim)
non_overlaps_idx <- which(str_detect(sequences, paste(rep("N", 10), collapse = "")))

# if there are any sequences that do not overlap,
# see if they can be collapsed and sum their
# counts

if (length(non_overlaps_idx) > 0) {

  non_overlaps <- sequences[non_overlaps_idx]
  non_overlaps_split = strsplit(non_overlaps, paste(rep("N", 10), collapse = ""))

  df_non_overlap <- data.frame(
    column_indexes = non_overlaps_idx,
    left_sequences = sapply(non_overlaps_split, "[[", 1),
    right_sequences = sapply(non_overlaps_split, "[[", 2)
  )

  ordered_left_sequences <- df_non_overlap %>%
    arrange(left_sequences, str_length(left_sequences)) %>%
    pull(left_sequences, name = column_indexes)

  idx <- 1
  # store modified left sequences in this list
  processed_left_sequences <- c()

  while (idx <= length(ordered_left_sequences)) {
    base_sequence <- ordered_left_sequences[idx]
    # 'x' is all sequences that match the base sequence. they should
    # already be together because of the prerequiste seqtab sorting above
    x <- ordered_left_sequences[
      startsWith(ordered_left_sequences, base_sequence)]

    t <- nchar(x)
    if (length(t) == 2 && (max(t) - min(t) == 1)) {
      x <- substr(x, 1, min(t))
    }

    processed_left_sequences <- c(processed_left_sequences, x)

    # update index
    idx <- idx + length(x)
  }

  # Now do the same for the right hand sequences!

  # reverse the sequences for easier comparison
  tmp <- data.frame(
    column_indexes = df_non_overlap$column_indexes,
    reversed = sapply(lapply(strsplit(df_non_overlap$right_sequences, NULL), rev), paste, collapse="")
  )
  ordered_reversed_right_sequences <- tmp %>%
    arrange(reversed, str_length(reversed)) %>%
    pull(reversed, name = column_indexes)

  idx <- 1
  processed_right_sequences <- c()

  while (idx <= length(ordered_reversed_right_sequences)) {
    base_sequence <- ordered_reversed_right_sequences[idx]

    x <- ordered_reversed_right_sequences[
      startsWith(ordered_reversed_right_sequences, base_sequence)]

    t <- nchar(x)
    if (length(t) == 2 && (max(t) - min(t) == 1)) {
      x <- substr(x, 1, min(t))
    }

    processed_right_sequences <- c(processed_right_sequences, x)

    # update index
    idx <- idx + length(x)
  }

  processed_right_sequences <- sapply(lapply(strsplit(processed_right_sequences, NULL), rev), paste, collapse="")

  # re-order the sequences by their indexes (ie. [1..n])
  processed_left_sequences <- processed_left_sequences[
    order(as.integer(names(processed_left_sequences)))]

  processed_right_sequences <- processed_right_sequences[
    order(as.integer(names(processed_right_sequences)))]

  # combine processed sequences
  df_non_overlap <- cbind(df_non_overlap,
                          processed_left_sequences = processed_left_sequences,
                          processed_right_sequences = processed_right_sequences
  )

  # combine the processed sequences
  combined <- NULL
  for (idx in 1:nrow(df_non_overlap)) {
    combined <- c(combined, paste(c(df_non_overlap[idx, ]$processed_left_sequences, df_non_overlap[idx, ]$processed_right_sequences), collapse = ""))
  }
  df_non_overlap$combined <- combined

  saveRDS(df_non_overlap,file="non_overlapping_seqs.RDS")
  write.table(df_non_overlap,file="non_overlapping_seqs.txt",quote=F,sep="\t",col.names=T,row.names=F)

  # finally write back to sequences
  sequences[df_non_overlap$column_indexes] <- df_non_overlap$combined
}

## Reference table for sequences - this is to reduce memory usage in intermediate files
df.sequences <- data.frame(
  sid = sprintf("S%d", 1:length(sequences)),
  sequences = sequences
)

## Setup parallel backend if asked

if (args$parallel) {
  n_cores <- ifelse(args$n_cores <= 0, detectCores(), args$n_cores)
  registerDoMC(n_cores)
} else {
  registerDoSEQ()
}


## II. Check for homopolymers.

if (!is.null(args$homopolymer_threshold) && args$homopolymer_threshold > 0) {

  ref_sequences <- readDNAStringSet(args$refseq_fasta)
  ref_names <- unique(names(ref_sequences))

  sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)

  # This object contains the aligned ASV sequences
  df_aln <- NULL
  df_aln <- foreach(seq1 = 1:length(sequences), .combine = "rbind") %dopar% {
    aln <- pairwiseAlignment(ref_sequences, sequences[seq1], substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE)
    num <- which.max(score(aln))
    patt <- c(alignedPattern(aln[num]), alignedSubject(aln[num]))
    ind <- sum(str_count(as.character(patt),"-"))
    data.frame(
      sid = df.sequences[seq1,]$sid,
      hapseq = as.character(patt)[2],
      refseq = as.character(patt)[1],
      refid = names(patt)[1],
      score = score(aln[num]),
      indels = ind
    )
  }

  # saveRDS(df_aln,file="alignments.RDS")
  # write.table(df_aln,file="alignments.txt",quote=F,sep="\t",col.names=T,row.names=F)

  df_aln <- df_aln %>% filter(score > args$alignment_threshold)

  masked_sequences <- readDNAStringSet(args$masked_fasta)
  df_masked <- NULL

  df_masked <- foreach(seq1 = 1:nrow(df_aln), .combine = "rbind") %dopar% {

    seq_1 <- DNAString(df_aln$refseq[seq1])

    asv_prime <-  DNAString(df_aln$hapseq[seq1])
    ref_rle <- Rle(as.vector(seq_1))

    # additionally, look for potential homopolymers that would exist in our
    # pairwise alignment reference sequence if it were not for insertions
    # in the sample haplotype in those regions.

    ref_rge <- ranges(ref_rle)
    mask_ranges <- IRanges()
    for (dna_base in DNA_ALPHABET[1:4]) {
      dna_ranges <- reduce(ref_rge[runValue(ref_rle) == dna_base | runValue(ref_rle) == "-"])

      vb <- NULL
      for (i in 1:length(dna_ranges)) {
        vb <- c(vb, sum(as.vector(seq_1[dna_ranges[i]]) %in% dna_base) > args$homopolymer_threshold)
      }

      mask_ranges <- append(mask_ranges, dna_ranges[vb])
    }

    maskseq <- getSeq(masked_sequences, df_aln$refid[seq1])
    trseq <- DNAString(as.character(maskseq))
    tr_rle <- Rle(as.vector(trseq))

    tr_rge <- ranges(tr_rle)[runValue(tr_rle) == "N"]
    gap_range <- ref_rge[runValue(ref_rle) == "-"]

    if (length(tr_rge) > 0) {
      for (i in 1:length(tr_rge)) {
        n_inserts <- sum(as.vector(seq_1[1:start(tr_rge[i])]) == '-')
        gaps_over_tr_rge <- gap_range[gap_range %over% shift(tr_rge[i], n_inserts)]
        n_over_range_gaps <- ifelse(length(gaps_over_tr_rge) == 0, 0, width(gaps_over_tr_rge))
        new_range <- IRanges(start = start(tr_rge[i]) + n_inserts, end = end(tr_rge[i]) + n_inserts + n_over_range_gaps)

        seq_1[new_range] <- "N" # mask refseq prime
        if (length(gaps_over_tr_rge) > 0) {
          for (j in 1:length(gaps_over_tr_rge)) {
            seq_1[gaps_over_tr_rge[j]] <- "-"
          }
        }
      }
    }

    tr_masks <- IRanges()
    ref_rle <- Rle(as.vector(seq_1))
    ref_rge <- ranges(ref_rle)
    tr_ranges <- reduce(ref_rge[runValue(ref_rle) == 'N' | runValue(ref_rle) == "-"])

    if (length(tr_ranges) > 0) {
      for (i in 1:length(tr_ranges)) {
        vseq <- as.vector(seq_1[tr_ranges[i]])
        if (all(c('N', '-') %in% vseq) || 'N' %in% vseq) {
          mask_ranges <- append(mask_ranges, tr_ranges[i])
        }
      }
    }

    mask_ranges <- unique(sort(mask_ranges))
    mask_ranges <- reduce(mask_ranges)
    reference_ranges <- mask_ranges

    pos <- c(DNA_ALPHABET[1:4], "N")
    if (length(mask_ranges) > 0) {
      for (i in 1:length(mask_ranges)) {

        vec_subseq <- as.vector(seq_1[reference_ranges[i]])
        ibase <- which(DNA_ALPHABET[1:4] %in% vec_subseq)
        dna_base <- ifelse(length(ibase) == 0, "N", DNA_ALPHABET[1:4][ibase])

        n_base <- sum(vec_subseq == dna_base)
        n_gaps <- sum(vec_subseq == '-')
        mask_ranges[i] <- mask_ranges[i] + sum(dna_base != 'N')

        replacement <- NULL

        # for homopolymers, replacement should be N * length of the mask
        # for tandem repeats, replacement should be length of original N mask
        trun <- end(mask_ranges[i]) - length(asv_prime)
        trun <- ifelse(trun < 0, 0, trun)
        if (dna_base %in% DNA_ALPHABET[1:4]) {
          replacement <- paste(as(Rle('N', width(mask_ranges[i]) - n_gaps - trun), 'character'), collapse = "")
        } else {
          replacement <- paste(as(Rle('N', n_base - trun), 'character'), collapse = "")
        }

        left_str <- ifelse(start(mask_ranges[i]) <= 1, "", as(subseq(asv_prime, 1, start(mask_ranges[i]) - 1), "character"))
        right_str <- ifelse(end(mask_ranges[i]) >= length(asv_prime), "", as(subseq(asv_prime, end(mask_ranges[i]) + 1, length(asv_prime)), "character"))
        asv_prime <- DNAString(paste(c(left_str, replacement, right_str), collapse = ""))

        mask_ranges <- shift(mask_ranges, -n_gaps)
      }
    }

    data.frame(
      sid = df_aln[seq1, ]$sid,
      refid = df_aln[seq1, ]$refid,
      asv_prime = as.character(asv_prime)
    )
  }

  # df_seqs <- inner_join(df_aln, df_masked, by = c("original", "refid", "refseq", "hapseq"))

  seqtab.nochim.df <- tibble::rownames_to_column(as.data.frame(t(seqtab.nochim)), "sequences") %>%
    inner_join(df.sequences, by = c("sequences")) %>%
    select(-c(sequences))

  df_seqs <- df_aln %>%
    ungroup() %>%
    select(sid, refid) %>%
    distinct() %>%
    inner_join(
      df_masked %>%
        select(sid, refid, asv_prime) %>%
        distinct()
      , by = c("sid", "refid"))

  seqtab.nochim.df <- df_seqs %>%
    inner_join(seqtab.nochim.df, by = c("sid")) %>%
    group_by(refid, asv_prime) %>%
    summarise(across(-c(sid), sum)) %>%
    ungroup() %>%
    select(-c(refid))

  seqtab.nochim.df <- as.data.frame(seqtab.nochim.df)
  rownames(seqtab.nochim.df) <- seqtab.nochim.df$asv_prime
  seqtab.nochim.df <- seqtab.nochim.df %>% select(-c(asv_prime))

  seqtab.nochim.df <- as.data.frame(t(seqtab.nochim.df))

  seqtab.nochim.df <- tibble::rownames_to_column(seqtab.nochim.df, "sample")
  seqtab.nochim.df <- seqtab.nochim.df %>% arrange(sample)
  seqtab.nochim.df[seqtab.nochim.df==0]=NA


} else {
  seqtab.nochim.df = as.data.frame(seqtab.nochim)
  seqtab.nochim.df$sample = rownames(seqtab.nochim)
  seqtab.nochim.df[seqtab.nochim.df==0]=NA
}

print("Done masking sequences...")

pat="-1A_|-1B_|-1_|-2_|-1AB_|-1B2_"
seqtab.nochim.df = seqtab.nochim.df %>%
  pivot_longer(cols = -c(sample), names_to = "asv",values_to = "reads",values_drop_na=TRUE) %>%
  mutate(locus = paste0(sapply(strsplit(sample,"_"),"[",1),"_",sapply(strsplit(sample,"_"),"[",2),"_",sapply(strsplit(sample,"_"),"[",3)))%>%
  mutate(sampleID = sapply(strsplit(sapply(strsplit(sample,pat),"[",2),"_trimmed"),"[",1)) %>%
  select(sampleID,locus,asv,reads)


temp = seqtab.nochim.df %>% select(locus,asv) %>% distinct()
loci =unique(temp$locus)
k=1
allele.sequences = data.frame(locus = seq(1,nrow(temp)),allele = seq(1,nrow(temp)),sequence = seq(1,nrow(temp)))
for(i in seq(1,length(loci))){
  temp2 = temp %>% filter(locus==loci[i])
  for(j in seq(1,nrow(temp2))){
    allele.sequences$locus[k+j-1] = loci[i]
    allele.sequences$allele[k+j-1] = paste0(loci[i],".",j)
    allele.sequences$sequence[k+j-1] = temp2$asv[j]
  }
  k=k+nrow(temp2)
}

allele.data = seqtab.nochim.df %>%
  left_join(allele.sequences %>% select(-locus),by=c("asv"="sequence")) %>%
  group_by(sampleID,locus,allele) %>%
  mutate(norm.reads.allele = reads/sum(reads))%>%
  group_by(sampleID,locus) %>%
  mutate(norm.reads.locus = reads/sum(reads))%>%
  mutate(n.alleles = n()) %>%
  arrange(sampleID, locus, reads)

head(allele.data)

saveRDS(allele.data,file="allele_data.RDS")
write.table(allele.data,file="allele_data.txt",quote=F,sep="\t",col.names=T,row.names=F)

# get memory footprint of environment
out <- as.data.frame(sort( sapply(ls(),function(x){object.size(get(x))})))
colnames(out) <- "bytes"
out <- out %>% dplyr::mutate(MB = bytes / 1e6, GB = bytes / 1e9)
write.csv(out, "postproc_memory_profile.csv")
