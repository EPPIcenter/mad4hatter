library(logger)
log_threshold(WARN)
log_appender(appender_console)

load_library <- function(library_name) {
  output <- capture.output({
    suppressWarnings({
      library(library_name, character.only = TRUE)
    })
  }, type = "message")
  
  # Separate warnings from messages
  warnings <- warnings()
  
  # Log messages
  if(length(output) > 0) {
    log_info(paste("Message from", library_name, ":", paste(output, collapse = "; ")))
  }
  
  # Log warnings
  if(length(warnings) > 0) {
    log_warn(paste("Warning from", library_name, ":", paste(warnings, collapse = "; ")))
  }
}

load_library("argparse")
load_library("dplyr")
load_library("magrittr")
load_library("stringr")

parser <- ArgumentParser(description='Handles denoising of sequences formed by joining reads through concatenation by DADA2')
parser$add_argument('--clusters', type="character", required=TRUE, help="RDS Clusters from DADA2. This is the main output from the DADA module.")
parser$add_argument('--log-level', type="character", default = "INFO", help = "Log level. Default is INFO.")

args <- parser$parse_args()

# Set the logging threshold 
log_level_arg <- match.arg(args$log_level, c("DEBUG", "INFO", "WARN", "ERROR", "FATAL"))
log_threshold(log_level_arg)

# Log the arguments parsed
args_string <- paste(sapply(names(args), function(name) {
  paste(name, ":", args[[name]])
}), collapse = ", ")

log_debug(paste("Arguments parsed successfully:", args_string))

## FOR DEBUGGING
# args=list()
# setwd("/home/bpalmer/Documents/GitHub/mad4hatter/work/e8/fbf5c2d9377a2e58bb236de5c872d7")
# args$clusters="dada2.clusters.txt"

# load the output from the dada2 process (allele_data saved into dada2.clusters.RDS)
clusters=NULL
if (grepl(".RDS|.rds", args$clusters)) {
  clusters=readRDS(args$clusters)
} else {
  clusters=read.table(args$clusters, header=T)
}

## I. Check for non overlapping sequences.

# make data frame clusters that has the distinct asvs and records in which row (sample/locus/allele in clusters) it showed
clusters %<>% dplyr::mutate(seqid=row_number())

clusters.1 = clusters %>%
  group_by(locus,asv) %>%
  dplyr::mutate(seqid=paste(seqid,collapse=";")) %>%
  select(locus,asv,seqid) %>%
  distinct()

# make a list of sequences
sequences = clusters.1$asv
# identify the sequences that didn't overlap as those that have 10 N's, which is how they are concatenated in merge in dada2
non_overlaps_idx <- which(str_detect(sequences, paste(rep("N", 10), collapse = "")))


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

  df_non_overlap %<>%
    mutate(combined = paste0(processed_left_sequences,"N",processed_right_sequences))

  write.table(df_non_overlap,file="concatenated_sequences.txt",quote=F,sep="\t",col.names=T,row.names=F)

  # finally write back to sequences
  sequences[df_non_overlap$column_indexes] <- df_non_overlap$combined
}

clusters.1$asv <- sequences
clusters.1 %<>%
  group_by(locus,asv) %>%
  dplyr::mutate(seqid=strsplit(seqid, ";")) %>%
  tidyr::unnest(cols=c(seqid)) %>%
  dplyr::mutate(seqid=as.integer(seqid))

clusters %<>%
    select(-c("asv")) %>%
    inner_join(clusters.1, by=c("locus","seqid")) %>%
    group_by(sampleID,locus,asv) %>%
    dplyr::mutate(reads=sum(reads)) %>%
    group_by(sampleID,locus) %>%
    mutate(norm.reads.locus = reads/sum(reads))%>%
    mutate(n.alleles = n()) %>%
    ungroup() %>%
	  distinct()

clusters$seqid=NULL
write.table(clusters,file="clusters.concatenated.collapsed.txt",quote=F,sep="\t",col.names=T,row.names=F)
