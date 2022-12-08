library(tidyverse)


args = commandArgs(trailingOnly=T)

alleledata_FILE=args[1]
codontable_FILE=args[2]
res_markers_info_FILE=args[3]
inputDIR=args[4]

##Allele table from dada2##

#alleledata_FILE="/home/isglobal.lan/ddatta/Projects/Pipeline/Runs/220909_M07977_0002_000000000-K7C6V/allele_data.txt"
alleledata=read.delim(alleledata_FILE,header=T,sep="\t")

##File containing amino acid code from DNA triplet##

#codontable_FILE="/home/isglobal.lan/ddatta/Projects/Pipeline/Files/codontable.txt"
codontable=read.delim(codontable_FILE,header=F,sep="\t")

##File containing the amplicon infos for the resistance marker positions##

#res_markers_info_FILE="/home/isglobal.lan/ddatta/Projects/Pipeline/Files/resistance_markers_amplicon_v4.txt"
res_markers_info=read.delim(res_markers_info_FILE,header=T,sep="\t")
res_markers_info=res_markers_info %>% filter(Codon_Start > 0, Codon_Start < ampInsert_length) %>% distinct(V5, .keep_all = T)

##Directory containing the mpileup of all the alleles##

#inputDIR="/home/isglobal.lan/ddatta/Projects/Pipeline/Runs/220909_M07977_0002_000000000-K7C6V/Mapping"
mpileupfiles<- list.files(path=inputDIR,pattern="*mpileup.txt$",full.names=T)
allelenames=str_replace_all(basename(mpileupfiles),".mpileup.txt","")
ampliconnames=sapply(strsplit(allelenames, "\\."), `[`, 1)

pat="-1A$|-1B$|-1AB$|-1B2$|-2$"

##Matrix with rows consisting of all the alleles, while columns consisting of the different resistance markers##
mat_refcodon=matrix("",nrow=length(allelenames),ncol=nrow(res_markers_info))
mat_altcodon=matrix("",nrow=length(allelenames),ncol=nrow(res_markers_info))

rownames(mat_refcodon)=allelenames
colnames(mat_refcodon)=res_markers_info$V5

rownames(mat_altcodon)=allelenames
colnames(mat_altcodon)=res_markers_info$V5

##Loop over each resistance marker##

for ( ii in 1:nrow(res_markers_info) )
{
resmarker=res_markers_info$V5[ii]
gene_strand=res_markers_info$V4[ii]
amplicon=sapply(strsplit(res_markers_info$amplicon[ii],pat),"[",1)
ampliconmpileupfiles=mpileupfiles[str_detect(mpileupfiles, amplicon, negate = FALSE)] ##allele mpileup files for the resistance marker##
print(amplicon)
print(ampliconmpileupfiles)
if(length(ampliconmpileupfiles) >0)
{
for ( jj in 1:length(ampliconmpileupfiles))  ##Calculate for each relevant allele##
{
allele_name=str_replace_all(basename(ampliconmpileupfiles[jj]),".mpileup.txt","")
temp1=read.delim(ampliconmpileupfiles[jj],header=F)
temp1=temp1 %>% mutate(V5=case_when(V4 =="<*>" ~ V3, TRUE ~ sapply(strsplit(V4, "\\,"), `[`, 1))) %>% data.frame()  ##If Reference, mpileup is a "*", else its Alternate,"*" ##

dna_triplet=temp1 %>% dplyr::filter(V2 %in% c(res_markers_info$Codon_Start[ii]:(res_markers_info$Codon_Start[ii]+2))) %>% select(V3,V5)  ##Obtain the dna triplet and codon change##
refdna_triplet=paste(dna_triplet$V3,collapse="")
if( gene_strand == "-" )  {refdna_triplet=intToUtf8(rev(utf8ToInt(chartr("ATGC", "TACG", refdna_triplet ))))}   ##Reverse complement the sequence if the gene is in the negative strand##
if(refdna_triplet!="")   ##Amplicon doesn't span the resistance marker##
{
refdna_codon=codontable[which(codontable[,1]==refdna_triplet),3]
mat_refcodon[allele_name,resmarker]=refdna_codon
}

altdna_triplet=paste(dna_triplet$V5,collapse="")
if( gene_strand == "-" )  {altdna_triplet=intToUtf8(rev(utf8ToInt(chartr("ATGC", "TACG", altdna_triplet ))))}   ##Reverse complement the sequence if the gene is in the negative strand##
if(altdna_triplet!="")   ##Amplicon doesn't span the resistance marker##
{
altdna_codon=codontable[which(codontable[,1]==altdna_triplet),3]
mat_altcodon[allele_name,resmarker]=altdna_codon
}

}
}
}


mat_refcodon=data.frame(mat_refcodon)
mat_altcodon=data.frame(mat_altcodon)


##################Haplotypes#################

 multi_amplicons= res_markers_info %>% dplyr::count(amplicon) %>% filter(n>1) %>% select(amplicon)   ##Alleles which span multiple amplicon##
 res_markers_multi_amplicons=res_markers_info %>% filter( amplicon %in% multi_amplicons$amplicon ) %>% data.frame()
 
 pat="-1A$|-1B$|-1AB$|-1B2$|-2$"
 
 ##Matrix with cols being the resistance markers which can form a haplotype and rows the alleles spanning them ( basically a subset of mat_altcodon )##
 altcodon_haps=mat_altcodon %>% filter(sapply(strsplit(rownames(mat_altcodon),"\\."),"[",1) %in% res_markers_multi_amplicons$amplicon) %>% select(str_replace_all(res_markers_multi_amplicons$V5,"-","."))

 arr_altcodon_haps=rep("",nrow(altcodon_haps))  ##Array with haplotypes for each allele which spans multiple amplicons##
 resmarker_haps=rep("",nrow(altcodon_haps))
 
 for(kk1 in 1:nrow(altcodon_haps))
 {
  codons_to_combine=res_markers_multi_amplicons %>% filter(amplicon %in% strsplit(rownames(altcodon_haps)[kk1],"\\.")[[1]]) %>% mutate(V5=str_replace_all(V5,"-",".")) %>% select(V5)
  resmarker_haps[kk1]=paste(str_replace_all(codons_to_combine$V5,"PF3D7_[0-9]+\\.",""),collapse="/")
  arr_altcodon_haps[kk1]=paste(altcodon_haps[kk1,codons_to_combine$V5],collapse="/")
 }
 
names(arr_altcodon_haps)=rownames(altcodon_haps)

mat_altcodon_haps= matrix("",nrow(altcodon_haps),length(unique(resmarker_haps)))   ##Matrix with cols being the haplotypes and rows the alleles spanning them##
mat_altcodon_haps=data.frame(mat_altcodon_haps)
rownames(mat_altcodon_haps)=rownames(altcodon_haps)
colnames(mat_altcodon_haps)=unique(resmarker_haps)

for(ii in 1:nrow(mat_altcodon_haps))
{
 mat_altcodon_haps[rownames(altcodon_haps)[ii],resmarker_haps[ii]]=arr_altcodon_haps[ii]
}
 
 ##Matrix with rows consisting of the samples, while columns consisting of the different resistance marker haplotypes##
 
sampleslist=unique(alleledata$sampleID)

 sample_altcodon_haps=matrix("",length(sampleslist),ncol(mat_altcodon_haps))
 sample_altcodon_haps=data.frame(sample_altcodon_haps)
 
 sample_altcodon_haps_reads=sample_altcodon_haps
 
##Loop over each sample##

for(dd1 in 1:length(sampleslist))
{
samplealleles=alleledata %>% filter(sampleID %in% sampleslist[dd1]) %>% select(allele,reads) %>% data.frame()  ##Alleles for the sample##

indmatch=which(rownames(mat_altcodon_haps) %in% samplealleles$allele)
if(length(indmatch) > 0)
{
sample_altcodon_haps[dd1,]=mat_altcodon_haps[indmatch,] %>%  dplyr::summarise(across(everything(), ~ paste(unique(.x[.x!=""]),collapse=","))) 

temp2_alt_haps_reads=samplealleles[rep(2, each = ncol(mat_altcodon_haps))]
colnames(temp2_alt_haps_reads)=colnames(mat_altcodon_haps)
rownames(temp2_alt_haps_reads)=samplealleles$allele
temp2_alt_haps_reads=temp2_alt_haps_reads %>% dplyr::filter(rownames(temp2_alt_haps_reads) %in% rownames(mat_altcodon_haps))
gg1=mat_altcodon_haps %>% dplyr::filter(rownames(mat_altcodon_haps) %in% samplealleles$allele) %>% mutate_all(~ case_when(. != "" ~ 1, TRUE ~ 0)) %>% data.frame()
temp2_alt_haps_reads2=temp2_alt_haps_reads*as.numeric(unlist(gg1)) 
temp2_alt_haps_reads2=temp2_alt_haps_reads2 %>%  dplyr::summarise(across(everything(), ~ paste(unique(.x[.x!=0]),collapse=","))) 
rownames(temp2_alt_haps_reads2)=sampleslist[dd1]
sample_altcodon_haps_reads[dd1,]=temp2_alt_haps_reads2

}
}

rownames(sample_altcodon_haps)=sampleslist
colnames(sample_altcodon_haps)=colnames(mat_altcodon_haps)
sample_altcodon_haps <- rownames_to_column(sample_altcodon_haps, "SampleName")

rownames(sample_altcodon_haps_reads)=sampleslist
colnames(sample_altcodon_haps_reads)=colnames(mat_altcodon_haps)
sample_altcodon_haps_reads <- rownames_to_column(sample_altcodon_haps_reads, "SampleName")

sample_altcodon_haps_FINAL=sample_altcodon_haps
sample_altcodon_haps_FINAL[-1] <- sprintf('%s (%s)', as.matrix(sample_altcodon_haps[-1]), as.matrix(sample_altcodon_haps_reads[-1]))
sample_altcodon_haps_FINAL=sample_altcodon_haps_FINAL %>%  mutate_all(funs(str_replace_all(., "\\(\\)", "")))

write.table(sample_altcodon_haps_FINAL,file="resmarkers_haplotype_summary.txt",quote=F,sep="\t",col.names=T,row.names=F)


##################Individual resistance markers#################

##Matrix with rows consisting of the samples, while columns consisting of the different resistance markers##

sampleslist=unique(alleledata$sampleID)

sample_refcodon=""
sample_altcodon=""

sample_altcodon_reads=""



##Loop over each sample##

for(kk in 1:length(sampleslist))
{
samplealleles=alleledata %>% filter(sampleID %in% sampleslist[kk]) %>% select(allele,reads) %>% data.frame()  ##Alleles for the sample##

##Ref Allele##
temp2_ref=mat_refcodon %>% dplyr::filter(rownames(mat_refcodon) %in% samplealleles$allele) %>%  dplyr::summarise(across(everything(), ~ paste(unique(.x[.x!=""]),collapse=","))) 
rownames(temp2_ref)=sampleslist[kk]
sample_refcodon=rbind(sample_refcodon,temp2_ref)

##Observed Allele##
temp2_alt=mat_altcodon %>% dplyr::filter(rownames(mat_altcodon) %in% samplealleles$allele) %>%  dplyr::summarise(across(everything(), ~ paste(unique(.x[.x!=""]),collapse=","))) 
rownames(temp2_alt)=sampleslist[kk]
sample_altcodon=rbind(sample_altcodon,temp2_alt)

temp2_alt_reads=samplealleles[rep(2, each = ncol(mat_altcodon))]
colnames(temp2_alt_reads)=colnames(mat_altcodon)
rownames(temp2_alt_reads)=samplealleles$allele
gg2=mat_altcodon %>% dplyr::filter(rownames(mat_altcodon) %in% samplealleles$allele) %>% mutate_all(~ case_when(. != "" ~ 1, TRUE ~ 0)) %>% data.frame()
temp2_alt_reads2=temp2_alt_reads*as.numeric(unlist(gg2)) 
temp2_alt_reads2=temp2_alt_reads2 %>%  dplyr::summarise(across(everything(), ~ paste(unique(.x[.x!=0]),collapse=","))) 
rownames(temp2_alt_reads2)=sampleslist[kk]
sample_altcodon_reads=rbind(sample_altcodon_reads,temp2_alt_reads2)
}

sample_refcodon=sample_refcodon[-1,]
colnames(sample_refcodon)=res_markers_info$V5
sample_refcodon <- rownames_to_column(sample_refcodon, "SampleName") %>% data.frame()

sample_altcodon=sample_altcodon[-1,]
colnames(sample_altcodon)=res_markers_info$V5
sample_altcodon <- rownames_to_column(sample_altcodon, "SampleName") %>% data.frame()

sample_altcodon_reads=sample_altcodon_reads[-1,]
colnames(sample_altcodon_reads)=res_markers_info$V5
sample_altcodon_reads <- rownames_to_column(sample_altcodon_reads, "SampleName") %>% data.frame()

sample_altcodon_FINAL=sample_altcodon
sample_altcodon_FINAL[-1] <- sprintf('%s (%s)', as.matrix(sample_altcodon[-1]), as.matrix(sample_altcodon_reads[-1]))
sample_altcodon_FINAL=sample_altcodon_FINAL %>%  mutate_all(funs(str_replace_all(., "\\(\\)", "")))

write.table(sample_altcodon_FINAL,file="resmarkers_summary.txt",quote=F,sep="\t",col.names=T,row.names=F)







