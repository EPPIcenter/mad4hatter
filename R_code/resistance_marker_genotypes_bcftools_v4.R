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

if (args$parallel) {
  n_cores <- ifelse(args$n_cores <= 0, detectCores(), args$n_cores)
  registerDoMC(n_cores)
} else {
  registerDoSEQ()
}



##Allele table from dada2##
allele.data=read.delim(args$alleledata_FILE,header=T,sep="\t")

##File containing amino acid code from DNA triplet##
codon.table=read.delim(args$codontable_FILE,header=F,sep="\t")

# load the reference sequence from fasta file
refseq_dnaset=readDNAStringSet(args$refseq)

# load the information for the resistance markers we want to take a look at
res_markers_info=read.delim(args$res_markers_info_FILE,header=T,sep="\t") %>%
  filter(Codon_Start > 0, Codon_Start < ampInsert_length) %>% distinct(V5, .keep_all = T) 

#simplify the allele table to keep the locus and the cigar string and only for those loci that are in the resistance marker table
unique_asvs = allele.data  %>% 
  ungroup() %>% 
  distinct(locus, pseudo_cigar_simple)  %>% 
  filter(locus %in% res_markers_info$amplicon)

refseqs = tibble(locus = names(refseq_dnaset),refseq = as.character(refseq_dnaset))  %>% 
  filter(locus %in% unique_asvs$locus & locus %in% res_markers_info$amplicon)

#simplify the table to keep the marker, its locus and the codon start, and to keep only the ones that have any data
res_markers_info_simple =res_markers_info%>%
  select(V5,amplicon,Codon_Start,V4) %>% 
  dplyr::rename(marker=V5,locus=amplicon,orientation=V4) %>% 
  filter(locus %in% unique_asvs$locus) 


# for each snp, get the allele table for that locus and join it with the resistance marker table
reversecomplement_pseudoCIGAR = function(string){
    pseudo_cigar_digits = strsplit(string, "\\D+")
    pseudo_cigar_letters = map(strsplit(string, "\\d+"), ~ .x[.x != ""])
    pseudo_cigar_split=Map(paste, pseudo_cigar_digits, pseudo_cigar_letters, sep = "")
    reversed = paste(rev(unlist(pseudo_cigar_split)),collapse="")
    complemented = chartr("TACG","ATGC",reversed)
  return(complemented)
}

res_markers_alleles = res_markers_info_simple %>% 
  left_join(refseqs, by="locus")   %>% 
  mutate(refseq_orientation = ifelse(orientation=="-", lapply(refseq, function(x) as.character(reverseComplement(DNAString(x)))), refseq),
         Codon_Start = ifelse(orientation=="-", nchar(refseq)-(Codon_Start + 2) +1,Codon_Start),
         reference_codon = substr(refseq_orientation, Codon_Start, Codon_Start+2)) %>% 
  left_join(codon.table  %>% select(V1,V2), by=c("reference_codon"="V1")) %>% 
  dplyr::rename(reference_aa=V2) %>%
  select(-refseq,-refseq_orientation) %>% 
  left_join(unique_asvs, by="locus")  %>% 
  mutate(pseudo_cigar_simple_rc = unlist(ifelse(orientation =="-",lapply(pseudo_cigar_simple, reversecomplement_pseudoCIGAR),pseudo_cigar_simple)) ) 
  # I'm here. I need to make sure reference codon is taken from the right one and that the mutations are taken from the reversed complemented pseudo cigar string
res_markers_alleles_allmatch = res_markers_alleles  %>% 
  filter(grepl("^\\d+M$",pseudo_cigar_simple_rc)) %>% 
  mutate(codon=reference_codon,
    aa = reference_aa,
    codon_refalt = "REF",
    aa_refalt = "REF")

modify_codon = function(codon_pseudocigar,codon){
  codon_split = strsplit(codon,"")[[1]]
  codon_pseudocigar_split = strsplit(codon_pseudocigar,"")[[1]]
  codon_split[codon_pseudocigar_split!="M"] = codon_pseudocigar_split[codon_pseudocigar_split!="M"]
  codon = paste(codon_split,collapse="")
  return(codon)
  }


res_markers_alleles_no_allmatch = res_markers_alleles  %>% 
  filter(!grepl("^\\d+M$",pseudo_cigar_simple_rc)) %>% 
  mutate(
    pseudo_cigar_noins = gsub("\\d+I", "", pseudo_cigar_simple_rc),
    pseudo_cigar_digits = strsplit(pseudo_cigar_noins, "\\D+"),
    pseudo_cigar_letters = map(strsplit(pseudo_cigar_noins, "\\d+"), ~ .x[.x != ""]),
    pseudo_cigar_digits_cumsum = lapply(strsplit(pseudo_cigar_noins, "\\D+"), function(y) cumsum(c(0, y))),   
    idx1 = mapply(function(x,y) max(which(y>x)), pseudo_cigar_digits_cumsum, Codon_Start),
    idx2 = mapply(function(x,y) max(which(y>x)), pseudo_cigar_digits_cumsum, Codon_Start+1),
    idx3 = mapply(function(x,y) max(which(y>x)), pseudo_cigar_digits_cumsum, Codon_Start+2),
    cod1 = mapply(function(x,y) unlist(x)[y],pseudo_cigar_letters,idx1),
    cod2 = mapply(function(x,y) unlist(x)[y],pseudo_cigar_letters,idx2),
    cod3 = mapply(function(x,y) unlist(x)[y],pseudo_cigar_letters,idx3),
    codon_pseudocigar = paste0(cod1,cod2,cod3),
    codon_refalt = ifelse(codon_pseudocigar=="MMM","REF","ALT"),
    codon = ifelse(codon_refalt=="REF",reference_codon,mapply(modify_codon,codon_pseudocigar,reference_codon))
    )    %>% 
    select(colnames(res_markers_alleles_allmatch %>% select(-aa,-aa_refalt))) 

res_markers_alleles_no_allmatch_nochange = res_markers_alleles_no_allmatch  %>% 
  filter(codon_refalt=="REF") %>% 
  mutate(aa = reference_aa,
         aa_refalt = "REF")

res_markers_alleles_no_allmatch_change = res_markers_alleles_no_allmatch  %>% 
  filter(codon_refalt=="ALT")  %>% 
  left_join(codon.table  %>% select(V1,V2), by=c("codon"="V1")) %>%
  dplyr::rename(aa=V2) %>%
  mutate(aa_refalt = ifelse(aa==reference_aa,"REF","ALT"))

res_markers_alleles_all = rbind(
  res_markers_alleles_allmatch,
  res_markers_alleles_no_allmatch_nochange,
  res_markers_alleles_no_allmatch_change
)

allele_data_snps_raw = allele.data %>% 
  left_join(res_markers_alleles_all %>% 
    select(locus,pseudo_cigar_simple,marker,reference_codon,codon,codon_refalt,reference_aa,aa,aa_refalt),
    by = c("locus","pseudo_cigar_simple")) %>% 
  filter(locus %in% res_markers_info$amplicon) %>% 
  mutate(geneID = sapply(strsplit(marker,"-"),"[",1),
    gene = sapply(strsplit(marker,"-"),"[",2),
    codonID = sapply(strsplit(marker,"-"),tail,1)) 
    
allele_data_snps = allele_data_snps_raw %>% 
  select(sampleID,geneID,gene,codonID,reference_codon,codon,codon_refalt,reference_aa,aa,aa_refalt,reads)


out_allele_data_snps = allele_data_snps
colnames(out_allele_data_snps) = 
  c("sampleID","Gene_ID","Gene","Codon_ID","Reference_Codon","Codon","Codon_Ref/Alt","Reference_AA","AA","AA_Ref/Alt","Reads")

allele_data_snps_collapsed = allele_data_snps  %>% 
  group_by(sampleID,geneID,gene,codonID,reference_codon,codon,codon_refalt,reference_aa,aa,aa_refalt) %>% 
  summarize(reads = sum(reads))

out_allele_data_snps_collapsed = allele_data_snps_collapsed
colnames(out_allele_data_snps_collapsed) = 
  c("sampleID","Gene_ID","Gene","Codon_ID","Reference_Codon","Codon","Codon_Ref/Alt","Reference_AA","AA","AA_Ref/Alt","Reads")

allele_data_microhap = allele_data_snps_raw %>% 
  select(sampleID,locus,pseudo_cigar_simple,geneID,gene,codonID,reference_aa,aa,reads)  %>% 
  group_by(sampleID,locus,pseudo_cigar_simple,geneID,gene,reads) %>% 
  summarize(microhap_idx = paste(codonID,collapse="/"),
            microhap = paste(aa,collapse="/"),
            microhap_ref = paste(reference_aa,collapse="/")) %>% 
  mutate(microhap_refalt = ifelse(microhap==microhap_ref,"REF","ALT"))              

out_allele_data_microhap = allele_data_microhap %>% 
  select(sampleID,geneID,gene,microhap_idx,microhap_ref,microhap,microhap_refalt,reads)
colnames(out_allele_data_microhap) = 
  c("sampleID","Gene_ID","Gene","Microhaplotype_Index","Reference_Microhaplotype","Microhaplotype","Microhaplotype_Ref/Alt","Reads")


allele_data_microhap_collapsed = allele_data_microhap  %>% 
  group_by(sampleID,geneID,gene,microhap_idx,microhap_ref,microhap,microhap_refalt) %>% 
  summarize(reads = sum(reads))

out_allele_data_microhap_collapsed = allele_data_microhap_collapsed %>% 
  select(sampleID,geneID,gene,microhap_idx,microhap_ref,microhap,microhap_refalt,reads)
colnames(out_allele_data_microhap_collapsed) =  
  c("sampleID","Gene_ID","Gene","Microhaplotype_Index","Reference_Microhaplotype","Microhaplotype_Ref/Alt","Reads")

##File containing the amplicon infos for the resistance marker positions##
write.table(out_allele_data_snps_collapsed, file="resmarker_table.txt", quote = F, row.names = F,sep="\t")
write.table(out_allele_data_microhap_collapsed, file="resmarker_microhap_table.txt", quote = F, row.names = F,sep="\t")



# need to account for snps that show up in more than 1 amplicon
# need to account for deletions and insertions
# need to put back the "novel" SNPs. And change the naming because they may not me novel. Maybe call otherSNP?