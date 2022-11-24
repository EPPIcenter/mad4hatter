library(tidyverse)
library(stringr)
library(dplyr)
library(dada2)
library(foreach)
library(parallel)
library(doMC)
library(muscle)
library(BSgenome)

args = commandArgs(trailingOnly=T)
numargs=length(args)
verbose=FALSE # todo: make that work
# masked_fasta=[numargs - 3]
homopolymer_threshold = as.integer(args[numargs - 2])
refseq_fasta=args[numargs - 1]
masked_fasta=args[numargs]
load (args[1])

## I. Check for non overlapping sequences. 

sequences = colnames(seqtab.nochim)
non_overlaps_idx <- which(str_detect(sequences, paste(rep("N", 10), collapse = "")))

# if there are any sequences that do not overlap,
# see if they can be collapsed and sum their
# counts

if (length(non_overlaps_idx) > 0) {

    # Calculate the difference in length between
    # real sequences, and see if the delta sequence
    # has no conflicting bases at each position. If 
    # they do, overwrite with N's and collapse the 
    # resolved sequences. Else if each position is unique,
    # assume that the sequences are part of the resolved
    # sequence. 

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
    base_sequence <- ordered_left_sequences[idx]

    # store modified left sequences in this list
    processed_left_sequences <- c()

    while (idx <= length(ordered_left_sequences)) {
        # 'x' is all sequences that match the base sequence. they should
        # already be together because of the prerequiste seqtab sorting above
        x <- ordered_left_sequences[
            startsWith(ordered_left_sequences, base_sequence)]

        # 't' is table of lengths of unique sequences.
        # if there are counts greater than 1 for unique sequences
        # larger than the base sequence, then there must be empty
        # spaces at those position. Thus, the actually base cannot
        # be resolved between the unique sequences. 
        t <- table(nchar(unique(x)))

        # if there is more than one unique sequence that is longer
        # than the base sequence, then replace with N's starting 
        # at the end position. This includes other sequences longer 
        # than these sequences. 
        if (any(t > t[1])) {
            deltas <- nchar(x) - nchar(base_sequence) + 1
            substr(x, deltas, nchar(x)) <- rep("N", deltas)
        } else { # there can only be one t[length(t)] here with above if
            # note: replace with substr()
            m <- which.max(nchar(x))
            n <- names(x) # cache names
            x <- rep(x[m], length(x))
            names(x) <- n
        }

        # update processed sequences list
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
    base_sequence <- ordered_reversed_right_sequences[1]
    processed_right_sequences <- NULL

    while (idx <= length(ordered_reversed_right_sequences)) {
        x <- ordered_reversed_right_sequences[
            startsWith(ordered_reversed_right_sequences, base_sequence)]

        t <- table(nchar(unique(x)))

        if (any(t > t[1])) {
            deltas <- nchar(x) - nchar(base_sequence) + 1
            substr(x, deltas, nchar(x)) <- rep("N", deltas)
        } else { # there can only be one t[length(t)] here with above if
            # note: replace with substr()
            m <- which.max(nchar(x))
            n <- names(x) # cache names
            x <- rep(x[m], length(x))
            names(x) <- n
        }

        # unreverse the strings
        x <- sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

        processed_right_sequences <- c(processed_right_sequences, x)
        
        # update index
        idx <- idx + length(x)
    }

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
    df_non_overlap <- unite(
        df_non_overlap,
        combined,
        c(
            processed_left_sequences,
            processed_right_sequences
        ), sep="")

    # finally write back to sequences
    sequences[df_non_overlap$column_indexes] <- df_non_overlap$combined
}


## II. Check for homopolymers. 

ref_sequences <- readDNAStringSet(refseq_fasta)
ref_names <- unique(names(ref_sequences))

sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)

registerDoMC(detectCores())
df_aln <- NULL
df_aln <- foreach(seq1 = 1:length(sequences), .combine = "rbind") %dopar% {
  seq_1 <- sequences[seq1]
  aln <- pairwiseAlignment(ref_sequences, seq_1, substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE)
  num <- which.max(score(aln))
  patt <- c(alignedPattern(aln[num]), alignedSubject(aln[num]))
  dist <- adist(as.character(patt)[1], as.character(patt)[2])
  ind <- sum(str_count(as.character(patt),"-"))
  data.frame(
    original = seq_1,
    hapseq = as.character(patt)[2],
    refseq = as.character(patt)[1],
    refid = names(patt)[1],
    score = score(aln[num]),
    indels = ind
  )
}

masked_sequences <- readDNAStringSet(masked_fasta)
df_masked <- NULL
df_masked <- foreach(seq1 = 1:nrow(df_aln), .combine = "rbind") %dopar% {
  # for (seq1 in 1:nrow(df_final)) {
  seq_1 <- DNAString(df_aln$refseq[seq1])
  maskseq <- getSeq(masked_sequences, df_aln$refid[seq1])

  ref_rle <- Rle(as.vector(seq_1))
  ss <- ranges(ref_rle)[runLength(ref_rle) > homopolymer_threshold]
  
  # cover the flanks
  start(ss) <- start(ss) - 1
  end(ss) <- end(ss) + 1
  
  # mask the 
  x <- DNAString(as.character(maskseq))
  x <- Rle(as.vector(x))
  ss <- append(ss, ranges(x)[runValue(x) == "N"])

  
  if (length(ss) == 0) {
    # just return the haplotype
    return (
      data.frame(
        refid = df_aln[seq1, ]$refid,
        refseq = df_aln[seq1, ]$refseq,
        hapseq = df_aln[seq1, ]$hapseq,
        asv_prime = df_aln[seq1, ]$hapseq
      )
    )
  }
  
  for (idx in 1:length(ss)) {
    css <- ss[idx]
    if (start(css) <= 0) {
      start(css) <- 1
    }
    if (end(css) > length(seq_1) ) {
      end(css) <- length(seq_1)
    }
    seq_1[start(css) : end(css)] <- "N"
  }
  
  data.frame(
    refid = df_aln$refid[seq1],
    refseq = df_aln[seq1, ]$refseq,
    hapseq = df_aln[seq1, ]$hapseq,
    asv_prime = as.character(seq_1)
  )
}

df_seqs <- inner_join(df_aln, df_masked, by = c("refid", "refseq", "hapseq"))
seqtab.nochim.df <- as.data.frame(t(seqtab.nochim))

seqtab.nochim.df$original <- base::rownames(seqtab.nochim.df)
df_seqs <- inner_join(df_seqs, seqtab.nochim.df, by = "original")

# everything here makes sense

df_final <- df_seqs %>%
  pivot_longer(
    cols = -c(original, hapseq, refseq, refid, score, indels, asv_prime),
    names_to = "sample",
    values_to = "counts"
  )

# todo: reorder the columns so that sample is first

seqtab.nochim.df <- df_final %>%
  group_by(sample, asv_prime) %>%
  summarise(counts = sum(counts)) %>%
  pivot_wider(names_from = asv_prime, values_from = counts) %>%
  # filter(counts != 0) %>%
  ungroup()

print(length(colnames(seqtab.nochim.df[, 2:ncol(seqtab.nochim.df)])))
print(colnames(seqtab.nochim.df[, 2:ncol(seqtab.nochim.df)]))

x <- colSums(seqtab.nochim.df[, 2:ncol(seqtab.nochim.df)]) != 0
seqtab.nochim.df <- seqtab.nochim.df[, c(TRUE, x)] # want the first 'sample' column (could write differently)
seqtab.nochim.df[seqtab.nochim.df==0]=NA

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
  mutate(n.alleles = n())

saveRDS(allele.data,file="allele_data.RDS")
write.table(allele.data,file="allele_data.txt",quote=F,sep="\t",col.names=T,row.names=F)
