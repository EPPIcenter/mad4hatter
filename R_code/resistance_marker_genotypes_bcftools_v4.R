library(tidyverse)
library(argparse)
library(Biostrings)
library(BSgenome)
library(doMC)
library(foreach)

parser <- ArgumentParser(description='MAD4HATTER Haplotype SNP Annotator')
parser$add_argument('--alleledata_FILE', type="character",
                    help='Allele data from MAD4HATTER pipeline.', required=TRUE)
parser$add_argument('--codontable_FILE', type="character", required=TRUE, help="Codon table reference. Example can be found and used in the MAD4HATTER repo 'templates/' folder.")
parser$add_argument('--res_markers_info_FILE', type="character", help="Table of resistance markers that are of interest",required=TRUE)
parser$add_argument('--refseq', type="character", help="Reference that was used during Post-Processing in the MAD4HATTER Pipeline.")
parser$add_argument('--parallel', action='store_true', help="Whether to run with multiple cores")
parser$add_argument('--n-cores', type = 'integer', default = -1, help = "Number of cores to use. Ignored if running parallel flag is unset.")

args <- parser$parse_args()
print(args)

## FOR DEBUGGING
# setwd("/home/bpalmer/Documents/GitHub/mad4hatter/work/a9/e820f3b8c34728b24003717619f9c3")
# args=list()
# args$alleledata_FILE="allele_data.txt"
# args$codontable_FILE="codontable.txt"
# args$res_markers_info_FILE="resistance_markers_amplicon_v4.txt"
# args$refseq="v4_refseq.fasta"

pat="-1A$|-1B$|-1AB$|-1B2$|-2$"
##Allele table from dada2##
allele.data=read.delim(args$alleledata_FILE,header=T,sep="\t")

##File containing amino acid code from DNA triplet##
codon.table=read.delim(args$codontable_FILE,header=F,sep="\t")

ss.refseq=readDNAStringSet(args$refseq)

## Create a table that combines the provided resmarker information and the allele table.
## They will be joined by locus (many-to-many) relationship. Each row will contain
## an allele from a sample at that locus, and will be inspected for codon changes
## described in the resmarker table.

res_markers_info=read.delim(args$res_markers_info_FILE,header=T,sep="\t") %>%
  dplyr::filter(Codon_Start > 0, Codon_Start < ampInsert_length) %>% distinct(V5, .keep_all = T) %>%
  dplyr::inner_join(allele.data, by=c("amplicon"="locus"), multiple = "all") %>%
  dplyr::mutate(cigar.keys=str_extract_all(cigar, "I|S|D|M")) %>%
  dplyr::mutate(cigar.values=str_split(cigar, "I|S|D|M")) %>%
  filter(str_detect(cigar,"I|S|D"))

registerDoSEQ()
resmarker.table=foreach(ii = 1:nrow(res_markers_info), .combine="bind_rows") %dopar% {
  resmarker=res_markers_info[ii,]
  cigar.keys=unlist(resmarker$cigar.keys)
  cigar.values=unlist(resmarker$cigar.values)
  cigar.values=as.integer(cigar.values[cigar.values!=""]) # the last entry is always an empty string from the str_split, should figure this out...
  cigar.rle=Rle(cigar.keys, cigar.values)
  cigar.rge=ranges(cigar.rle)
  cigar.rge.1=cigar.rge[runValue(cigar.rle) != "D"]
  asv=DNAString(resmarker$asv)
  asv.1=asv[cigar.rge.1]

  codon.start=resmarker$Codon_Start
  strand=resmarker$V4
  codon=asv.1[codon.start:(codon.start+2)]

  if(strand=="-")  {
    codon=Biostrings::reverseComplement(codon)
  }

  ## Check for AA change
  codon=as.character(codon)
  resmarker$aa.expected=resmarker$V6
  resmarker$aa.actual=codon.table %>% filter(V1==codon) %>% pull(V2)
  resmarker$aa.change=resmarker$aa.actual!=resmarker$V6

  ## May not be necessary - should be optional. Check for codon change.
  ref=as.character(getSeq(ss.refseq,resmarker$amplicon))
  ref=DNAString(ref)
  ref.codon=ref[codon.start:(codon.start+2)]
  if(strand=="-")  {
    ref.codon=Biostrings::reverseComplement(ref.codon)
  }
  ref.codon=as.character(ref.codon)
  resmarker$codon.expected=ref.codon
  resmarker$codon.actual=codon
  resmarker$codon.change=codon != ref.codon

  return(resmarker)
}

## `V5` is the gene
resmarker.final=resmarker.table %>%
  select(sampleID, V5, Codon_Start, aa.actual, reads) %>%
  dplyr::rename(`locus_gene`=V5, `snp(AA)` = aa.actual, codon_start_position=Codon_Start)

##File containing the amplicon infos for the resistance marker positions##
write.table(resmarker.final, file="resmarker_table.txt", quote = F, row.names = F)



## Now look at novel SNPs
novel.snps=foreach(ii = 1:nrow(res_markers_info), .combine="bind_rows") %dopar% {
  resmarker=res_markers_info[ii,]
  asv=DNAString(resmarker$asv)

  ref=as.character(getSeq(ss.refseq,resmarker$amplicon))
  ref=DNAString(ref)

  cigar.keys=unlist(resmarker$cigar.keys)
  cigar.values=unlist(resmarker$cigar.values)
  cigar.values=as.integer(cigar.values[cigar.values!=""]) # the last entry is always an empty string from the str_split, should figure this out...
  cigar.rle=Rle(cigar.keys, cigar.values)
  cigar.rge=ranges(cigar.rle)

  cigar.rge.2=cigar.rge[runValue(cigar.rle) %in% c("S", "I", "D")]
  novel.cigar=NULL
  if (length(cigar.rge.2) > 0) {
    novel.cigar=foreach (jj = 1:length(cigar.rge.2), .combine="bind_rows") %do% {
      pos=paste(start(cigar.rge.2[jj]):end(cigar.rge.2[jj]), collapse=",")
      cigar=paste(rep(runValue(cigar.rle[cigar.rge.2[jj]]), width(cigar.rge.2[jj])), collapse = ",")

      tibble(pos, cigar)
    }

    novel.cigar=novel.cigar %>%
      dplyr::mutate(
        sampleID=resmarker$sampleID,
        gene_locus=resmarker$V5,
        codon_start=resmarker$Codon_Start,
        pos=strsplit(pos,","),
        cigar=strsplit(cigar,","),
        reads=resmarker$reads
      ) %>%
      tidyr::unnest(cols=c(pos,cigar)) %>%
      dplyr::mutate(
        pos=as.integer(pos),
        codon_start=as.integer(codon_start)
      )


    # only keep positions that are novel snps
    novel.cigar=novel.cigar %>%
      filter(!((pos >= codon_start) & (pos <= (codon_start+2) )))

    novel.cigar=foreach (jj = 1:nrow(novel.cigar), .combine="bind_rows") %do% {
      novel.cigar.row = novel.cigar[jj,]
      pos=as.integer(novel.cigar.row$pos)
      novel.cigar.row$REF=as.character(ref[pos])
      novel.cigar.row$ALT=as.character(asv[pos])

      novel.cigar.row
    }
  }
 return(novel.cigar)
}

novel.snps = novel.snps %>%
  select(pos, cigar, REF, ALT, reads, gene_locus, sampleID)

##File containing the amplicon infos for the resistance marker positions##
write.table(novel.snps, file="resmarker_novel_snps.txt", quote = F, row.names = F)
