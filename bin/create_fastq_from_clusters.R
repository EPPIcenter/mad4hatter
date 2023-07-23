library(argparse)

parser <- ArgumentParser(prog="Utility to convert DADA2 clusters object into fastqs by sample", description='DADA2 workflow script')
parser$add_argument('--clusters', type="character", help="Use quality scores from clusters to generate fastq files")
parser$add_argument('--output', type="character", help="Directory to write fastq files to")
parser$add_argument('--amplicon-table', type="character", help="Amplicon table - used to parse loci.")
args <- parser$parse_args()
print(args)

library(dplyr)
library(BSgenome)
library(stringr)

amplicon.table = args$amplicon_table


# create regex to extract sampleID (to be used with remove)
pat=paste(sprintf("^%s_", amplicon.table$amplicon), collapse="|") # find amplicons specified in amplicon table (with _)
pat=paste(c(pat, "_trimmed_merged.fastq.gz$"),collapse="|")  # remove _trimmed from end


## Write out the clusters
foreach::registerDoSEQ()
clusterx=foreach(ii=1:length(clusters), .combine="bind_rows") %do% {
  clusterID=names(clusters)[ii]
  clx=clusters[[clusterID]]$clustering
  locus=amplicon.table$amplicon[!is.na(str_extract(clusterID, amplicon.table$amplicon))]
  sampleID=str_remove_all(clusterID, pat)

  qs=as.integer(round(clusters[[clusterID]]$quality))
  qs=qs[!is.na(qs)]

  tbl=tibble(
    sampleID=sampleID,
    locus=locus,
    reads=clx$abundance,
    asv=clx$sequence, # identical to "denoised"
    quality=intToUtf8(qs + 33)
  )
}

clusterx=clusterx %>%
  group_by(sampleID, locus) %>%
  dplyr::mutate(n.alleles=n()) %>%
  inner_join(clusterx, by=c("sampleID", "locus", "reads", "asv", "quality")) %>%
  dplyr::mutate(allele=sprintf("%s.%d", locus, seq(n.alleles))) %>%
  select(sampleID, locus, allele, reads, asv, quality) %>%
  distinct()

write.table(clusterx,file="clusters.txt",quote = F, row.names = F)

if (!dir.exists(args$output)) {
  dir.create(args$output)
}

for (sid in unique(clusterx$sampleID)) {
  allele.data.filtered=allele.data[clusterx$sampleID==sid,]

  clusterID=sprintf("%s_trimmed_merged.fastq.gz", sid)
  clusters.filtered=clusters[[clusterID]]

  ss=DNAStringSet(allele.data.filtered$asv)
  names(ss)<-allele.data.filtered$allele

  qs = intToUtf8(clusters.filtered$quality) - 33


  qs=PhredQuality()

  out=QualityScaledDNAStringSet(ss,qs)
  writeQualityScaledXStringSet(out, filepath = file.path("cluster-fq", sprintf("%s.fastq", sid)))
}