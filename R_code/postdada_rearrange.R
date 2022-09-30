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
ampliconFILE=args[numargs]
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

df_final <- NULL
registerDoMC(detectCores())
df_final <- foreach(refidx = 1:length(ref_names), .combine = "rbind") %dopar% {

  ref_name <- ref_names[refidx]
  ref_seq <- getSeq(ref_sequences, ref_name)

  if (sum(str_detect(rownames(seqtab.nochim), ref_name)) == 0) {
    print(paste("WARN:", ref_name, "not found in sequence table!"))
    return(NULL)
  }

  df_subset <- seqtab.nochim[str_detect(rownames(seqtab.nochim), ref_name), ]
  if ((nrow(df_subset) == 1 && df_subset[1, colidx] == 0) || rowSums(df_subset) == 0) {
    print(paste("WARN:", ref_name, "could not be linked to any sequences!"))
    return(NULL)
  }

  colidx = unname(which(colSums(df_subset) > 0))
  sequences = colnames(df_subset)[colidx]

  set <- DNAStringSet(c(ref_seq, sequences))
  aln <- muscle(set, quiet = TRUE)

  ref_dna <- as.matrix(unmasked(aln)[1])
  ref_rle <- Rle(ref_dna)
  ref_rge <- ranges(ref_rle)


  for (dna_base in DNA_ALPHABET[1:4]) {

      # 1. find situations where we have 'AA--AA-A-A--'.
      # A homopolymer exists if the same base is found repeatedly around no base calls,
      # and the detected base count is greater than the set threshold. In the above 
      # example, A equals 7, thus the homopolymer extends through the entire string.


      indexes <- which(runValue(ref_rle) == "-" | runValue(ref_rle) == dna_base)
      irange <- IRanges(indexes)
      
      # find ranges where "-" and dna_base are next to each other
      sequential <- irange[width(irange) > 1]
      if (length(sequential) > 0) {
          for (i in seq.int(1,length(sequential))) {
              sub_rge <- ref_rge[start(sequential[i]):end(sequential[i])]
              sub_rge_start <- start(sub_rge[1])
              sub_rge_end <- end(sub_rge[length(sub_rge)])

              sub_rge_dna <- ref_dna[sub_rge_start:sub_rge_end]                
              homopolymer_len <- sum(sub_rge_dna == dna_base)
              if (homopolymer_len > homopolymer_threshold) {

                  # exclude possible variation next to homopolymers
                  if (sub_rge_start > 1)
                      sub_rge_start <- sub_rge_start + 1

                  if (sub_rge_end < length(ref_dna))
                      sub_rge_end <- sub_rge_end - 1

                  ref_dna[sub_rge_start:sub_rge_end] <- "N" # "N" is homopolymer mask
              }
          }
      }

      # 2. find repeating base counts (ie. 'AAAAA')
      base_run <- ref_rge[runValue(ref_rle) == dna_base]
      base_run_target_indexes <- which(width(base_run) > homopolymer_threshold)

      if (length(base_run_target_indexes) > 0) {
          for (idx in base_run_target_indexes) {

              # exclude possible variation next to homopolymers
              base_run_start <- start(base_run[idx])
              base_run_end <- end(base_run[idx])
              lower_idx <- ifelse(base_run_start > 1, base_run_start - 1, base_run_start)
              upper_idx <- ifelse(base_run_end < length(ref_dna), base_run_end + 1, base_run_end)
              ref_dna[lower_idx:upper_idx] <- "N" # "N" is homopolymer mask 

          }
      }
  } # end loop


  u_aln <- unmasked(aln)
  writeXStringSet(u_aln, paste0(ref_name, ".fasta"), append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

  asv_set <- NULL 
  for (idx in 2:nrow(aln)) {
      seq_mat <- as.matrix(u_aln[idx])
      seq_mat[which(ref_dna == "N")] <- "N" 
      seq_mat[which(ref_dna == "-")] <- "N" # mask no call
      seq_str <- paste(as.character(seq_mat), collapse="")
      seq_str <- str_replace_all(seq_str, "-", "N")
      asv_set <- c(asv_set, seq_str)
  }

  df_amplicon <- NULL
  for (cidx in 1:length(colidx)) {
    col <- colidx[cidx]
    seq <- asv_set[cidx]
    for (ridx in 1:nrow(df_subset)) {
        df_amplicon <- rbind(df_amplicon,
            data.frame(
              # amplicon = ref_name,
              sample = rownames(df_subset)[ridx],
              ref_sequences = seq,
              counts = unname(df_subset[ridx, col])
            ))
    }
  }
  return(filter(df_amplicon, counts != 0))
}

save(df_final, file="~/df_final.rda")

## III. Collapse sequences 

seqtab.nochim.df <- df_final %>%
    group_by(sample, ref_sequences) %>%
    summarise(counts = sum(counts)) %>%
    pivot_wider(names_from = ref_sequences, values_from = counts) %>%
    ungroup()

## IV. Write allele tables

# seqtab.nochim.df = as.data.frame(seqtab.nochim)
# seqtab.nochim.df$sample = rownames(seqtab.nochim)
# seqtab.nochim.df[seqtab.nochim.df==0]=NA
pat="-1A_|-1B_|-1_|-2_|-1AB_|-1B2_"
seqtab.nochim.df = seqtab.nochim.df %>% 
  pivot_longer(cols = -sample, names_to = "asv",values_to = "reads",values_drop_na=TRUE) %>% 
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