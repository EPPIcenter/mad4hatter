library(argparse)

parser <- ArgumentParser(description='Post process dada2 inferred sequences')
parser$add_argument('--homopolymer-threshold', type="integer",
                   help='homopolymer threshold to begin masking')
parser$add_argument('--refseq-fasta', type="character")
parser$add_argument('--masked-fasta', type="character")
parser$add_argument('--alignment-threshold', type="integer", default = 60)
parser$add_argument('--parallel', action='store_true')
parser$add_argument('--n-cores', type = 'integer', default = -1, help = "Number of cores to use. Ignored if running parallel flag is unset.")
parser$add_argument('--sample-coverage', type="character", help = "Sample coverage file from QC to append sample coverage statistics that are P. falciparum specific.")
parser$add_argument('--amplicon-coverage', type="character", help = "Amplicon coverage file from QC to append amplicon coverage statistics that are P. falciparum specific.")
parser$add_argument('--amplicon-table', type="character", required=TRUE, help = "Amplicon table with primer pools. This is used to organize the sequence table by amplicon.")
parser$add_argument('--clusters', type="character", required=TRUE, help="RDS Clusters from DADA2. This is the main output from the DADA module.")

args <- parser$parse_args()
print(args)

library(stringr)
library(dplyr)
library(dada2)
library(foreach)
library(parallel)
library(BSgenome)
library(tidyr)
library(doMC)
library(tibble)
library(ggplot2)
library(Biostrings)
library(magrittr)

# setwd("/home/bpalmer/Documents/GitHub/mad4hatter/work/13/444d5f387c683912ce0855f452fce4")
# # # FOR DEBUGGING
# args <- list()
# args$homopolymer_threshold <- 5
# args$refseq_fasta <- "v4_refseq.fasta"
# args$masked_fasta <- "v4_refseq.fasta.2.7.7.80.10.25.3.mask"
# args$dada2_output <- "seqtab.nochim.RDS"
# args$parallel <- FALSE
# args$alignment_threshold <- 60
# args$clusters = "dada2.clusters.RDS"
# args$amplicon_table="v4_amplicon_info.tsv"
# # args$pseudo_fastq_output="pseudo-fastqs"


# load the output from the dada2 process (allele_data saved into dada2.clusters.RDS)
clusters=NULL
if (grepl(".RDS|.rds", args$clusters)) {
  clusters=readRDS(args$clusters)
} else {
  clusters=read.table(args$clusters, header=T)
}

clusters %<>% 
  mutate(sampleID = str_remove_all(sampleID, pat = "_trimmed"))

## I. Check for non overlapping sequences.

# make data frame clusters.1 that has the distinct asvs and records in which row (sample/locus/allele in clusters) it showed
clusters.1=clusters %>%
  ungroup() %>%
  dplyr::mutate(seqid=row_number()) %>%
  group_by(locus,asv) %>%
  dplyr::mutate(seqid=paste(seqid,collapse=";")) %>%
  select(locus,asv,seqid) %>%
  distinct() 
# make a list of sequences
sequences = clusters.1$asv
# identify the sequences that didn't overlap as those that have 10 N's, which is how they are concatenated in merge in dada2
non_overlaps_idx <- which(str_detect(sequences, paste(rep("N", 10), collapse = "")))

# read in the amplicon info tsv which has the name and genomic location
amplicon.table=read.table(args$amplicon_table,header=T)

# if there are any sequences that do not overlap,
# see if they can be collapsed and sum their
# counts

# here we will collapse all the non overlapping sequences that only differed in the length of the reads on either side
# the thinking is that those differences are driven by differences in quality of reads rather than actual variation


# Note for Brian: in df_aln we would have only unique sequences but there's no collapsing of the truncated sequences so those are repeated
if (length(non_overlaps_idx) > 0) {

  # extract the non overlapping sequences
  non_overlaps <- sequences[non_overlaps_idx]
  # split the two ends of the sequence
  non_overlaps_split = strsplit(non_overlaps, paste(rep("N", 10), collapse = ""))

  # put both ends of the sequence and the positions in the sequences list in a data frame
  df_non_overlap <- data.frame(
    column_indexes = non_overlaps_idx,  # note for Brian: why are these called column_indexes if they are from rows?
    left_sequences = sapply(non_overlaps_split, "[[", 1),
    right_sequences = sapply(non_overlaps_split, "[[", 2)
  )

  # extract and order the left sequences and order them alphabetically (and keep the indexes in sequence as names)
  ordered_left_sequences <- df_non_overlap %>%
    arrange(left_sequences, str_length(left_sequences)) %>%
    pull(left_sequences, name = column_indexes)

  idx <- 1
  # store modified left sequences in this list
  processed_left_sequences <- c()

  # for each sequence find any others that start with the same sequence 
  while (idx <= length(ordered_left_sequences)) {
    base_sequence <- ordered_left_sequences[idx]
    # 'x' is all sequences that match the base sequence. they should
    # already be together because of the prerequiste seqtab sorting above 
    # sequences were sorted alphabetically so matching ones will be together, and then sorted by length so the shortest will come first
    # note that some left sequences will be identical because the differences in the asv are in the right sequence
    x <- ordered_left_sequences[
      startsWith(ordered_left_sequences, base_sequence)]
    
    if (length(x) > 1) {     # Note for Brian: why only if you have 2 and if the difference is 1? I'm temporarily changing to more than 1 and at most 3, but that's an arbitrary number to not get rid of true trs
    # if there are sequences that are longer than the first by more than 3 bases then those are not analyzed/modified
    #if (length(t) == 2 && (max(t) - min(t) == 1)) {  
      t <- nchar(x)
      idx.close = which((t-t[1])<4)
      x <- substr(x[idx.close], 1, t[1])
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
  # extract and order the reversed right sequences and order them alphabetically (and keep the indexes in sequence as names)
  ordered_reversed_right_sequences <- tmp %>%
    arrange(reversed, str_length(reversed)) %>%
    pull(reversed, name = column_indexes)
  # start idx to loop over sequences
  idx <- 1
  # store modified sequences in this list
  processed_right_sequences <- c()

  while (idx <= length(ordered_reversed_right_sequences)) {
    base_sequence <- ordered_reversed_right_sequences[idx]

    x <- ordered_reversed_right_sequences[
      startsWith(ordered_reversed_right_sequences, base_sequence)]

    # if there are more than 1 sequences , take the ones that differ by less than 4 bases and shorten to the length of the first one, if there are longer ones, leave out for the next round
    if (length(x) >1) {
      t <- nchar(x)
      idx.close = which((t-t[1])<4)
      x <- substr(x[idx.close], 1, t[1])
    }

    #append the sequences, modified if so
    processed_right_sequences <- c(processed_right_sequences, x)

    # update index
    idx <- idx + length(x)
  }

  # reverse the sequences
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
  #combined <- NULL
  #for (idx in 1:nrow(df_non_overlap)) {
  #  combined <- c(combined, paste(c(df_non_overlap[idx, ]$processed_left_sequences, df_non_overlap[idx, ]$processed_right_sequences), collapse = ""))
  #}
  #df_non_overlap$combined <- combined

  df_non_overlap %<>% 
    mutate(combined = paste0(processed_left_sequences,"N",processed_right_sequences)) 

  saveRDS(df_non_overlap,file="non_overlapping_seqs.RDS")
  write.table(df_non_overlap,file="non_overlapping_seqs.txt",quote=F,sep="\t",col.names=T,row.names=F)

  # finally write back to sequences
  sequences[df_non_overlap$column_indexes] <- df_non_overlap$combined
}
clusters.1$sequences <- sequences
## Reference table for sequences - this is to reduce memory usage in intermediate files
## This table keeps track of which sample/locus had the unique sequence
df.sequences.idx <- data.frame(
  seqid = sprintf("S%d", 1:length(sequences)),
  asvid = clusters.1$seqid
)

## Setup parallel backend if asked

if (args$parallel) {
  n_cores <- ifelse(args$n_cores <= 0, detectCores(), args$n_cores)
  registerDoMC(n_cores)
} else {
  registerDoSEQ()
}

# read the sequences from the reference genome (already extracted into a fasta file)
ref_sequences <- readDNAStringSet(args$refseq_fasta)
# read the names (loci)
#ref_names <- unique(names(ref_sequences))  # note to Brian: commenting out as it doesn't show up again


# Note to Brian: please explain what this does and why those numbers are chosen
sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)

#########  IMPORTANT: NON PF SEQUENCES NEED TO BE INCLUDED IN REFERENCE!


####################### NEED TO ANNOTATE BETTER BELOW THIS   ##############################


# This object contains the aligned ASV sequences
df_aln <- NULL
df_aln <- foreach(seq1 = 1:nrow(clusters.1), .combine = "bind_rows") %dopar% {
  # the alignment is performed only with the reference sequence for the corresponding locus as we have this information from the demultiplexing step
  refseq.seq1 = ref_sequences[clusters.1$locus[seq1]] #commenting out next lines as I concatenated all genomes
  #if(substr(clusters.1$locus[seq1],1,2)=="Pf"){   # temporary fix to deal with other species
  #  refseq.seq1 = ref_sequences[clusters.1$locus[seq1]]
  #}else{
  #  refseq.seq1 = ref_sequences["Pf3D7_13_v3-1041593-1041860-1AB"]  # in the meantime I'll use Pf's ldh as reference, which is wrong because it's not the same sequence...!
  #}
  aln <- pairwiseAlignment(refseq.seq1, str_remove_all(sequences[seq1],"N"), substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE)
  patt <- c(alignedPattern(aln), alignedSubject(aln))
  ind <- sum(str_count(as.character(patt),"-"))
  data.frame(
    seqid = df.sequences.idx[seq1,]$seqid,
    hapseq = as.character(patt)[2],
    refseq = as.character(patt)[1],
    refid = clusters.1$locus[seq1],
    score = score(aln),
    indels = ind
  )
}

saveRDS(df_aln,file="alignments.RDS")
write.table(df_aln,file="alignments.txt",quote=F,sep="\t",col.names=T,row.names=F)

## add a histogram of the scores
## add a vline where the filter threshold is
# pdf(filename="alignments.pdf")
g = ggplot(df_aln) +
  geom_histogram(aes(df_aln$score)) +
  geom_vline(xintercept = args$alignment_threshold) +
  ggtitle("Distribution of alignment scores and the alignment threshold") +
  xlab("Alignment Score") +
  ylab("Frequency")

ggsave(filename="alignments.pdf", g, device="pdf")

## Filter off targets here
df_aln %<>% filter(score > args$alignment_threshold)

## Note for Brian: commenting the next 3 lines out because they don't show up again
# create regex to extract sampleID (to be used with `str_remove_all()`)
#pat=paste(sprintf("^%s_", amplicon.table$amplicon), collapse="|") # find amplicons specified in amplicon table (with _)
#pat=paste(c(pat, "_trimmed_merged.fastq.gz$"),collapse="|")  # remove _trimmed from end

## II. Check for homopolymers.

# This only needs to be done in the reference sequences. We may want to change how homopolymers are defined in the future. As in we could also use homopolymers observed in data, not in reference to define one that doesn't make the cut in the reference
 
if (!is.null(args$homopolymer_threshold) && args$homopolymer_threshold > 0) {
  # Define a function to process each sequence
  homomask_sequence <- function(sequence) {
    # somewhat convoluted way to get the Rle from the sequence in the DNAStringSet object that will be passed to the function
    ref_seq_rle <- Rle(as.vector(unlist(strsplit(as.character(sequence), ""))))
    # replace the homopolymers with Ns
    ref_seq_rle@values[ref_seq_rle@lengths > args$homopolymer_threshold] <- "N"
    # add 2 N's to the run to mask 1 base upstream and 1 base downstream
    ref_seq_rle@lengths[ref_seq_rle@lengths > args$homopolymer_threshold] <-  as.integer(ref_seq_rle@lengths[ref_seq_rle@lengths > args$homopolymer_threshold] + 2)
    # remove one N from the run to mask 1 base upstream and 1 base downstream
    ref_seq_rle@lengths[lag(ref_seq_rle@lengths > args$homopolymer_threshold,default=FALSE)] <-  as.integer(ref_seq_rle@lengths[lag(ref_seq_rle@lengths > args$homopolymer_threshold,default=FALSE)] -1)
    ref_seq_rle@lengths[lead(ref_seq_rle@lengths > args$homopolymer_threshold,default=FALSE)] <-  as.integer(ref_seq_rle@lengths[lead(ref_seq_rle@lengths > args$homopolymer_threshold,default=FALSE)] -1)
    # if a base was between homopolymers we have a -1 now, so we need to account by taking an N out of one of the homopolymers, and I'll do from the left
    ref_seq_rle@lengths[lead(ref_seq_rle@lengths == -1 ,default=FALSE)] <-  as.integer(ref_seq_rle@lengths[lead(ref_seq_rle@lengths == -1 ,default=FALSE)]  -1)
    # if homopolymer was at the beginning or end we added an N, so let's remove it
    ref_seq_rle@lengths[c(1, length(ref_seq_rle@lengths))] <- ref_seq_rle@lengths[c(1, length(ref_seq_rle@lengths))] - (ref_seq_rle@values[c(1, length(ref_seq_rle@lengths))] == "N")
    # remove empty runs, make -1's, 0
    ref_seq_rle <- Rle(ref_seq_rle@values, ifelse(ref_seq_rle@lengths<0, 0, ref_seq_rle@lengths))
    # paste the Rle back together into a sequence
    homomasked_seq <- DNAString(paste(rep(runValue(ref_seq_rle), times = runLength(ref_seq_rle)), collapse = ""))
    # return the sequence
    return(homomasked_seq)
  }
  # Apply the function to each element in ref_sequences
  homomask_ref_sequences <-DNAStringSet(lapply(ref_sequences, homomask_sequence))
}

# load the masked sequences:
trmask_ref_sequences = readDNAStringSet(args$masked_fasta)

# now we have ref_sequences, homomask_ref_sequences and trmask_ref_sequences. they all have the same number of objects and each object has the same length across them 

#Note to Brian:
#for now I'll keep the choice using args$homopolymer_threshold, but we should add the option to do tr masking only etc. And the homopolymers should come out of this code. 
#also if no masking is required skip the processes to mask. 

merge_homo_tr = function(seq.homo, seq.tr) {
  seq1=as.character(seq.homo)
  seq2=as.character(seq.tr)
  combined <- character(max(nchar(seq1), nchar(seq2)))
  combined[(strsplit(seq1,"")[[1]]=="N") | (strsplit(seq2,"")[[1]]=="N")] <- "N"
  combined[!((strsplit(seq1,"")[[1]]=="N") | (strsplit(seq2,"")[[1]]=="N"))] <- strsplit(seq1[[1]],"")[[1]][!((strsplit(seq1,"")[[1]]=="N") | (strsplit(seq2,"")[[1]]=="N"))]
  combined_rle = Rle(combined)
  combined_seq = DNAString(paste(rep(runValue(combined_rle), times = runLength(combined_rle)), collapse = ""))
  return(combined_seq)
}

homo_tr_mask_ref_sequences = DNAStringSet(mapply(merge_homo_tr,homomask_ref_sequences,trmask_ref_sequences))

df_homo_tr_mask_ref_sequences = tibble(homo_tr_ref_seq = as.character(homo_tr_mask_ref_sequences),refid = names(homo_tr_mask_ref_sequences))

# Now let's mask the aligned reference sequences 
# first make a data frame with all the unique aligned sequences and join with the data frame that has the masked references
df_aln_references = df_aln  %>% 
  select(refid,refseq) %>% 
  distinct() %>% 
  left_join(df_homo_tr_mask_ref_sequences,by="refid")


mask_aligned_sequence = function(homo_tr_ref_seq,refseq){
  string_homo_tr_ref_seq = homo_tr_ref_seq
  split_string <- gregexpr("[^-]+|-+", refseq)
  split_string <- regmatches(refseq, split_string)[[1]]
  newstring = character()
  for(i in split_string){
    if(grepl("-",i)){
      newstring=paste0(newstring,i)
    }else{
      newstring=paste0(newstring,substr(string_homo_tr_ref_seq,1,nchar(i)))
      string_homo_tr_ref_seq = substr(string_homo_tr_ref_seq,nchar(i)+1,nchar(string_homo_tr_ref_seq))
    }
  }
  return(newstring)
}

remove_dashes_between_N = function(newstring){
  NdashN = unique(regmatches(newstring,gregexpr("N-+N", newstring))[[1]])
  newstring_nodash = newstring
  for(i in NdashN){
    newstring_nodash=gsub(i,strrep("N",nchar(i)),newstring_nodash)
  }
  return(newstring_nodash)
  }



df_aln_references %<>% 
  mutate(masked_aligned_refseq = mapply(mask_aligned_sequence,homo_tr_ref_seq,refseq)) %>% 
  mutate(masked_aligned_nodashinN_refseq = sapply(masked_aligned_refseq,remove_dashes_between_N))


#update df_aln 
df_aln %<>% left_join(df_aln_references,by=c("refid","refseq"))


# now let's make a pseudoCIGAR

if (!is.null(args$homopolymer_threshold) && args$homopolymer_threshold > 0) {
  allele.data <- foreach(seq2cigar = 1:nrow(df_aln), .combine = "rbind") %dopar% {

      hapseq = df_aln[seq2cigar,]$hapseq
      masked_aligned_refseq=df_aln[seq2cigar,]$masked_aligned_refseq[[1]]

      bs.cigar=BString(x=paste(rep("M", nchar(hapseq)), collapse=""))

      hapseq_dna = DNAString(hapseq)
      hapseq_rle = Rle(as.vector(hapseq_dna))
      hapseq_rge = ranges(hapseq_rle)

      masked_aligned_refseq_dna = DNAString(masked_aligned_refseq)
      masked_aligned_refseq_rle = Rle(as.vector(masked_aligned_refseq_dna))
      masked_aligned_refseq_rge = ranges(masked_aligned_refseq_rle)

      masked_ranges = masked_aligned_refseq_rge[runValue(masked_aligned_refseq_rle) == "N"]

      # ranges_pseudocigar is where we'll add each of the ranges for the pseudocigar
      ranges_pseudocigar=list()

      # note that insertions are ambiguous!

      insertion_rge = masked_aligned_refseq_rge[runValue(masked_aligned_refseq_rle) == "-"]
      ranges_pseudocigar[["I"]]=insertion_rge[insertion_rge %outside% masked_ranges]

      deletion_rge = hapseq_rge[runValue(hapseq_rle) == "-"]
      ranges_pseudocigar[["D"]]<-deletion_rge[deletion_rge %outside% masked_ranges]

      ranges_pseudocigar[["N"]]<-masked_ranges

          for (dna_base in c("A", "C", "T", "G")) {
            refbase = masked_aligned_refseq_rge[runValue(masked_aligned_refseq_rle) == dna_base]
            hapbase = hapseq_rge[runValue(hapseq_rle) == dna_base]
            snprge = IRanges::setdiff(hapbase, refbase)
            snprge = snprge[snprge %outside% masked_ranges & snprge %outside% insertion_rge]
            ranges_pseudocigar[[dna_base]] <- snprge
          }


          for (ii in names(ranges_pseudocigar)) {
            if (length(ranges_pseudocigar[[ii]]) > 0) {
              for (jj in 1:length(ranges_pseudocigar[[ii]])) {
                bs.cigar[ranges_pseudocigar[[ii]][jj]] <- c(ii)
              }
            } else {
              bs.cigar[ranges_pseudocigar[[ii]]] <- c(ii)
            }
          }

          remove_insertions_between_N = function(newstring){
            NIN = unique(regmatches(newstring,gregexpr("NI+N", newstring))[[1]])
            newstring_noI = newstring
            for(i in NIN){
              newstring_noI=gsub(i,"NN",newstring_noI)
            }
            return(newstring_noI)
          }

          pseudo_cigar_build = remove_insertions_between_N(as.character(bs.cigar))

          pseudo_cigar_rle = Rle(strsplit(pseudo_cigar_build,"")[[1]])

          pseudo_cigar=paste(sprintf("%d%s", runLength(pseudo_cigar_rle), runValue(pseudo_cigar_rle)), collapse = "")

          pseudo_cigar_build_simple = str_replace_all(pseudo_cigar_build,"N","M")

          pseudo_cigar_simple_rle = Rle(strsplit(pseudo_cigar_build_simple,"")[[1]])

          pseudo_cigar_simple=paste(sprintf("%d%s", runLength(pseudo_cigar_simple_rle), runValue(pseudo_cigar_simple_rle)), collapse = "")


          data.frame(
            seqid = df_aln[seq2cigar, ]$seqid,
            refid = df_aln[seq2cigar, ]$refid,
            asv_prime = hapseq,
            pseudo_cigar_simple=pseudo_cigar_simple,
            pseudo_cigar = pseudo_cigar
          )
  }
}else{
  allele.data <- foreach(seq2cigar = 1:nrow(df_aln), .combine = "rbind") %dopar% {

      hapseq = df_aln[seq2cigar,]$hapseq
      refseq=df_aln[seq2cigar,]$refseq[[1]]

      bs.cigar=BString(x=paste(rep("M", nchar(hapseq)), collapse=""))

      hapseq_dna = DNAString(hapseq)
      hapseq_rle = Rle(as.vector(hapseq_dna))
      hapseq_rge = ranges(hapseq_rle)

      refseq_dna = DNAString(refseq)
      refseq_rle = Rle(as.vector(refseq_dna))
      refseq_rge = ranges(refseq_rle)

      # ranges_pseudocigar is where we'll add each of the ranges for the pseudocigar
      ranges_pseudocigar=list()

      # note that insertions are ambiguous!

      insertion_rge = refseq_rge[runValue(refseq_rle) == "-"]
      ranges_pseudocigar[["I"]]=insertion_rge

      deletion_rge = hapseq_rge[runValue(hapseq_rle) == "-"]
      ranges_pseudocigar[["D"]]<-deletion_rge


          for (dna_base in c("A", "C", "T", "G")) {
            refbase = refseq_rge[runValue(refseq_rle) == dna_base]
            hapbase = hapseq_rge[runValue(hapseq_rle) == dna_base]
            snprge = IRanges::setdiff(hapbase, refbase)
            snprge = snprge[snprge %outside% insertion_rge]
            ranges_pseudocigar[[dna_base]] <- snprge
          }


          for (ii in names(ranges_pseudocigar)) {
            if (length(ranges_pseudocigar[[ii]]) > 0) {
              for (jj in 1:length(ranges_pseudocigar[[ii]])) {
                bs.cigar[ranges_pseudocigar[[ii]][jj]] <- c(ii)
              }
            } else {
              bs.cigar[ranges_pseudocigar[[ii]]] <- c(ii)
            }
          }


          pseudo_cigar_build = as.character(bs.cigar)

          pseudo_cigar_rle = Rle(strsplit(pseudo_cigar_build,"")[[1]])

          pseudo_cigar=paste(sprintf("%d%s", runLength(pseudo_cigar_rle), runValue(pseudo_cigar_rle)), collapse = "")

          data.frame(
            seqid = df_aln[seq2cigar, ]$seqid,
            refid = df_aln[seq2cigar, ]$refid,
            asv_prime = hapseq,
            pseudo_cigar_simple=pseudo_cigar,
            pseudo_cigar = pseudo_cigar
          )
  }
}


  ## ASV in the original user file should be replaced with the identfier to use less memory
  allele.data %<>%
    dplyr::rename(locus=refid) %>%
    inner_join(df.sequences.idx, by=c("seqid")) %>%
    inner_join(clusters.1,by=c('locus','asvid'='seqid')) %>%
    inner_join(clusters,by=c('locus', "asv")) %>%
    select(sampleID, locus, asv_prime, reads, allele, pseudo_cigar_simple, pseudo_cigar) %>%
    group_by(sampleID, locus, asv_prime) %>%
    mutate(reads=sum(reads)) %>%
    dplyr::rename(asv=asv_prime) %>%
    distinct()

print("Done with pseudoCIGAR...")
print(allele.data[1,])
saveRDS(allele.data,file="allele_data.RDS")
write.table(allele.data,file="allele_data.txt",quote=F,sep="\t",col.names=T,row.names=F)

## QC Postprocessing

if (!is.null(args$sample_coverage) && file.exists(args$sample_coverage)) {
  sample.coverage <- read.table(args$sample_coverage, header = TRUE, sep = "\t") %>%
    pivot_wider(names_from = "X", values_from = "NumReads")

  qc.postproc <- sample.coverage %>% 
    left_join(clusters  %>% 
      ungroup()  %>% 
      select(sampleID,reads) %>% 
      group_by(sampleID) %>%
      summarise(OutputDada2 = sum(reads)), by = c("SampleName" = "sampleID")
    )  %>% 
    left_join(allele.data  %>% 
      ungroup()  %>% 
      select(sampleID,reads) %>% 
      group_by(sampleID) %>%
      summarise(OutputPostprocessing = sum(reads)), by = c("SampleName" = "sampleID")
    )  %>%
    pivot_longer(cols = c(Input, `No Dimers`, Amplicons, OutputDada2, OutputPostprocessing))

  colnames(qc.postproc) <- c("SampleName","","NumReads")
  write.table(qc.postproc, quote=F,sep='\t',col.names = TRUE, row.names = F, file = args$sample_coverage)
}

if (!is.null(args$amplicon_coverage) && file.exists(args$amplicon_coverage)) {
  amplicon.coverage <- read.table(args$amplicon_coverage, header = TRUE, sep = "\t")

  qc.postproc <- amplicon.coverage  %>% 
    left_join(clusters %>% 
      group_by(sampleID,locus) %>%
      summarise(OutputDada2 = sum(reads)),
      by = c("SampleName" = "sampleID", "Amplicon" = "locus"),
      ) %>% 
    left_join(allele.data %>%
      group_by(sampleID,locus) %>%
      summarise(OutputPostprocessing = sum(reads)),
          by = c("SampleName" = "sampleID", "Amplicon" = "locus")) 
   qc.postproc$OutputDada2[is.na(qc.postproc$OutputDada2)] <- 0
   qc.postproc$OutputPostprocessing[is.na(qc.postproc$OutputPostprocessing)] <- 0 
   
  write.table(qc.postproc, quote=F,sep='\t',col.names = TRUE, row.names = F, file = args$amplicon_coverage)
}

## get memory footprint of environment

out <- as.data.frame(sort( sapply(ls(),function(x){object.size(get(x))})))
colnames(out) <- "bytes"
out <- out %>% dplyr::mutate(MB = bytes / 1e6, GB = bytes / 1e9)
write.csv(out, "postproc_memory_profile.csv")
