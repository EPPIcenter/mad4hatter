library(tidyverse)


args = commandArgs(trailingOnly=T)

alleledata_FILE=args[1]
codontable_FILE=args[2]
res_markers_info_FILE=args[3]
inputDIR=args[4]

##Allele table from dada2##

#alleledata_FILE="/home/isglobal.lan/ddatta/Projects/Pipeline/ampseq_workflow_dd/22_0085_Test/allele_data.txt"
alleledata=read.delim(alleledata_FILE,header=T,sep="\t")

##File containing amino acid code from DNA triplet##

#codontable_FILE="/home/isglobal.lan/ddatta/Projects/Pipeline/ampseq_workflow_dd/codontable.txt"
codontable=read.delim(codontable_FILE,header=F,sep="\t")

##File containing the amplicon infos for the resistance marker positions##

#res_markers_info_FILE="/home/isglobal.lan/ddatta/Projects/Pipeline/ampseq_workflow_dd/resistance_markers_amplicon_v3.txt"
res_markers_info=read.delim(res_markers_info_FILE,header=T,sep="\t")
res_markers_info=res_markers_info %>% filter(Codon_Start > 0, Codon_Start < ampInsert_length) %>% distinct(V5, .keep_all = T)

##Directory containing the mpileup of all the alleles##

#inputDIR="/home/isglobal.lan/ddatta/Projects/Pipeline/ampseq_workflow_dd/22_0085_Test/Mapping3"
mpileupfiles<- list.files(path=inputDIR,pattern="*mpileup.txt$",full.names=T)
allelenames=str_replace_all(basename(mpileupfiles),".mpileup.txt","")
ampliconnames=sapply(strsplit(allelenames, "\\."), `[`, 1)

pat="-1A$|-1B$|-1AB$|-2$"

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
for ( jj in 1:length(ampliconmpileupfiles))  ##Calculate for each relevant allele##
{
allele_name=str_replace_all(basename(ampliconmpileupfiles[jj]),".mpileup.txt","")
temp1=read.delim(ampliconmpileupfiles[jj],header=F)
temp1=temp1 %>% mutate(V5=case_when(V4 =="<*>" ~ V3, TRUE ~ sapply(strsplit(V4, "\\,"), `[`, 1))) %>% data.frame()  ##If Reference, mpileup is a "*", else its Alternate,"*" ##

dna_triplet=temp1 %>% dplyr::filter(V2 %in% c(res_markers_info$Codon_Start[ii]:(res_markers_info$Codon_Start[ii]+2))) %>% select(V3,V5)  ##Obtain the dna triplet and codon change##
refdna_triplet=paste(dna_triplet$V3,collapse="")
if( gene_strand == "-" )  {refdna_triplet=intToUtf8(rev(utf8ToInt(chartr("ATGC", "TACG", refdna_triplet ))))}   ##Reverse complement the sequence if the gene is in the negative strand##
refdna_codon=codontable[which(codontable[,1]==refdna_triplet),3]
mat_refcodon[allele_name,resmarker]=refdna_codon

altdna_triplet=paste(dna_triplet$V5,collapse="")
if( gene_strand == "-" )  {altdna_triplet=intToUtf8(rev(utf8ToInt(chartr("ATGC", "TACG", altdna_triplet ))))}   ##Reverse complement the sequence if the gene is in the negative strand##
altdna_codon=codontable[which(codontable[,1]==altdna_triplet),3]
mat_altcodon[allele_name,resmarker]=altdna_codon

}
}


##Matrix with rows consisting of the samples, while columns consisting of the different resistance markers##
mat_refcodon=data.frame(mat_refcodon)
mat_altcodon=data.frame(mat_altcodon)

sampleslist=unique(alleledata$sampleID)

sample_refcodon=""
sample_altcodon=""


##Loop over each sample##

for(kk in 1:length(sampleslist))
{
samplealleles=alleledata %>% filter(sampleID %in% sampleslist[kk]) %>% select(allele) %>% data.frame()  ##Alleles for the sample##

##Ref Allele##
temp2_ref=mat_refcodon %>% dplyr::filter(rownames(mat_refcodon) %in% samplealleles$allele) %>%  dplyr::summarise(across(everything(), ~ paste(unique(.x[.x!=""]),collapse=","))) 
rownames(temp2_ref)=sampleslist[kk]
sample_refcodon=rbind(sample_refcodon,temp2_ref)

##Observed Allele##
temp2_alt=mat_altcodon %>% dplyr::filter(rownames(mat_altcodon) %in% samplealleles$allele) %>%  dplyr::summarise(across(everything(), ~ paste(unique(.x[.x!=""]),collapse=","))) 
rownames(temp2_alt)=sampleslist[kk]
sample_altcodon=rbind(sample_altcodon,temp2_alt)
}

sample_refcodon=sample_refcodon[-1,]
colnames(sample_refcodon)=res_markers_info$V5
sample_refcodon <- rownames_to_column(sample_refcodon, "SampleName") %>% data.frame()

sample_altcodon=sample_altcodon[-1,]
colnames(sample_altcodon)=res_markers_info$V5
sample_altcodon <- rownames_to_column(sample_altcodon, "SampleName") %>% data.frame()

 write.table(sample_altcodon,file="resmarkers_summary.txt",quote=F,sep="\t",col.names=T,row.names=F)




