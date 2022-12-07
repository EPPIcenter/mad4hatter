library(stringr)
library(dplyr)
library(dada2)
library(foreach)
library(parallel)
library(doMC)
library(muscle)
library(BSgenome)
library(argparse)
library(tidyr)

parser <- ArgumentParser(description='Post process dada2 inferred sequences')
parser$add_argument('--homopolymer-threshold', type="integer",
                   help='homopolymer threshold to begin masking')
parser$add_argument('--refseq-fasta', type="character")
parser$add_argument('--masked-fasta', type="character")
parser$add_argument('--dada2-output', type="character", required = TRUE)
# parser$add_argument('--parallel', type='boolean', default=FALSE)

args <- parser$parse_args()
load (args$dada2_output)

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

if (!is.null(args$homopolymer_threshold) && args$homopolymer_threshold > 0) {
  
  ref_sequences <- readDNAStringSet(args$refseq_fasta)
  ref_names <- unique(names(ref_sequences))
  
  sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
  
  registerDoMC(detectCores())
  df_aln <- NULL
  df_aln <- foreach(seq1 = 1:length(sequences), .combine = "rbind") %dopar% {
    seq_1 <- sequences[seq1]
    aln <- pairwiseAlignment(ref_sequences, seq_1, substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE)
    num <- which.max(score(aln))
    patt <- c(alignedPattern(aln[num]), alignedSubject(aln[num]))
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

  df_aln <- df_aln %>% filter(score > 0)
  
  masked_sequences <- readDNAStringSet(args$masked_fasta)
  df_masked <- NULL
  df_masked <- foreach(seq1 = 1:nrow(df_aln), .combine = "rbind") %dopar% {
    # for (seq1 in 1:nrow(df_aln)) {
    
    mask_ranges <- NULL
    
    # mask homopolymers, and the bases flanking them
    # input: reference sequence from pairwise alignment
    seq_1 <- DNAString(df_aln$refseq[seq1])
    asv_prime <-  DNAString(df_aln$hapseq[seq1])
    
    ref_rle <- Rle(as.vector(seq_1))
    
    # additionally, look for potential homopolymers that would exist in our 
    # pairwise alignment reference sequence if it were not for insertions 
    # in the sample haplotype in those regions. 
    ref_rge <- ranges(ref_rle)
    for (dna_base in DNA_ALPHABET[1:4]) {

      mask_ranges <- append(mask_ranges, ranges(ref_rle)[runValue(ref_rle) == dna_base & runLength(ref_rle) > args$homopolymer_threshold])

      indexes <- which(runValue(ref_rle) == "-" | runValue(ref_rle) == dna_base)
      irange <- IRanges(indexes)
      
      # find ranges where "-" and dna_base are next to each other
      # 'AAAA--A' = insertions in sample within homopolymer.
      # this is what we want to catch
      sequential <- irange[width(irange) > 2]  # eg. '-A-' or 'A-A'
      
      # move on if we don't have an 'A-' patterns
      if (length(sequential) == 0)
        next
      
      for (i in 1:length(sequential)) {
        # 'sub_rge' = sub-range
        # these are the indexes that contain dna base of interest and 
        # '-' that are next to each other in the original sequence
        sub_rge <- ref_rge[sequential[i]]
        sub_rge_dna <- seq_1[sub_rge]
        n_base <- sum(as.vector(sub_rge_dna) %in% dna_base)
        if (n_base > args$homopolymer_threshold) {
          # indels in the masked region should result in a mask that is the same
          # length as the homopolymer region in the reference
          # note: there is probably a better way to do this....
          replacement <- paste(as(Rle(dna_base, n_base), 'character'), collapse = "")
          left_str <- subseq(asv_prime, 1, start(sub_rge[1]) - 1)
          right_str <- subseq(asv_prime, end(sub_rge[length(sub_rge)]) + 1)
          asv_prime <- DNAString(paste(c(as(left_str, "character"), replacement, as(right_str, "character")), collapse = ""))
          
          # we masked above but this will let us mask the flanks using the logic below
          new_range <- IRanges(start = start(sub_rge)[1], end = end(sub_rge)[length(sub_rge)] - (length(sub_rge) - n_base))
          mask_ranges <- append(mask_ranges, new_range)
        }
      }
    }
    
    
    # as mentioned, add additional masking around homopolymers.
    # if out of range, this will be taken care of during the actual masking.
    if (!is.null(mask_ranges))
      start(mask_ranges) <- start(mask_ranges) - 1
    if (!is.null(mask_ranges))
      end(mask_ranges) <- end(mask_ranges) + 1
    
    
    # mask tandem repeats
    # input: masked sequences from the pipeline
    
    maskseq <- getSeq(masked_sequences, df_aln$refid[seq1])
    trseq <- DNAString(as.character(maskseq))
    tr_rle <- Rle(as.vector(trseq))
    tr_rge <- ranges(tr_rle)[runValue(tr_rle) == "N"]
    
    if (length(tr_rge) > 0) {
      for (i in 1:length(tr_rge)) {
        # count how many gaps there are in the reference sequence from insertions,
        # and adjust the tandem repeat masks accordingly
        ref_dna <- DNAString(df_aln$refseq[seq1])
        n_gaps <- str_count(as.character(ref_dna[1:start(tr_rge[i])]), "-")
        new_range <- IRanges(start = start(tr_rge[i]) + n_gaps, end = end(tr_rge[i]) + n_gaps)
        replacement <- paste(as(Rle('N', width(tr_rge[i])), 'character'), collapse = "")
        
        # adjust the mask like before so that the masked range 
        left_str <- subseq(asv_prime, 1, start(tr_rge[i]) - 1)
        right_str <- subseq(asv_prime, end(tr_rge[i]) + 1)
        asv_prime <- DNAString(paste(c(as(left_str, "character"), replacement, as(right_str, "character")), collapse = ""))
        
        # note: don't worry about appending the mask here because we technically did already,
        # and we don't need to worry about masking the sides of the mask
        # mask_ranges <- append(mask_ranges, new_range)
      }
    }
    
    # mask_ranges <- append(mask_ranges, ranges(trseq)[runValue(trseq) == "N"])
    
    # if there is nothing to mask, just return
    if (length(mask_ranges) == 0) {
      return (
        data.frame(
          original = df_aln[seq1, ]$original,
          refid = df_aln[seq1, ]$refid,
          refseq = df_aln[seq1, ]$refseq,
          hapseq = df_aln[seq1, ]$hapseq,
          asv_prime = as.character(asv_prime)
        )
      )
    }
    
    for (idx in 1:length(mask_ranges)) {
      css <- mask_ranges[idx]
      if (start(css) <= 0) {
        start(css) <- 1
      }
      if (end(css) > length(asv_prime) ) {
        end(css) <- length(asv_prime)
      }
      asv_prime[start(css) : end(css)] <- "N"
    }
    
    data.frame(
      original = df_aln[seq1, ]$original,
      refid = df_aln$refid[seq1],
      refseq = df_aln[seq1, ]$refseq,
      hapseq = df_aln[seq1, ]$hapseq,
      asv_prime = as.character(asv_prime)
    )
  }
  
  df_seqs <- inner_join(df_aln, df_masked, by = c("original", "refid", "refseq", "hapseq"))
    
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
  
  seqtab.nochim.df[seqtab.nochim.df==0]=NA
  
} else {
  seqtab.nochim.df = as.data.frame(seqtab.nochim)
  seqtab.nochim.df$sample = rownames(seqtab.nochim)
  seqtab.nochim.df[seqtab.nochim.df==0]=NA    
}


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