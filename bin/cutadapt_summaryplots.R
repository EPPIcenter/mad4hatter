library(tidyverse)
library(gridExtra)
library(ggbeeswarm)
library(rmarkdown)
library(knitr)

args = commandArgs(trailingOnly=T)

summaryFILE=args[1]
samplestatFILE=args[2]
ampliconFILE=args[3]
outDIR=args[4]

##FOR DEBUGGING
# setwd("/home/bpalmer/Documents/GitHub/mad4hatter/work/1f/0ca34bcad8061f443b16cc0aefa1ea")
# summaryFILE="amplicon_coverage.txt"
# samplestatFILE="sample_coverage.txt"
# ampliconFILE="v4_amplicon_info.tsv"
# outDIR="quality_report"

df=read.table(summaryFILE,header=T)
df$SampleID=as.factor(df$SampleID)
df$Locus=as.factor(df$Locus)
df = df %>% 
  mutate(SampleNumber=sapply(str_split(SampleID,'_S'),tail,1)) %>% 
  mutate(SampleID=sapply(str_split(SampleID,'_S(\\d+)'),head,1)) %>% 
  mutate(Pool=sapply(str_split(Locus,'-'),tail,1)) %>% 
  arrange(SampleID) %>% 
  data.frame()

samples = df %>% select(SampleID,SampleNumber) %>% distinct() %>% 
  arrange(SampleNumber)

df$SampleID=factor(df$SampleID,levels=unique(samples$SampleID))

amplicon_stats=df %>% select(-Pool) %>% pivot_wider(names_from = Locus, values_from = Reads) %>% data.frame()
write.table(amplicon_stats, file=paste(outDIR,"/amplicon_stats.txt",sep=""), quote=F, sep ="\t", col.names=T, row.names=F)

sample_amplicon_stats=df %>% group_by(SampleID,Pool) %>% dplyr::summarise(medianReads=median(Reads)) %>% pivot_wider(names_from = Pool, values_from = medianReads) %>% data.frame()

# Declare a variable to contain the pools detected. Each
# panel has a different number of pools, so make this dynamic.
pool_columns <- sprintf("Pool_%s", sort(unique(df$Pool)))
colnames(sample_amplicon_stats)=c("SampleID",pool_columns)

loci_stats = df %>% group_by(SampleID) %>% group_by(SampleID,Pool) %>% dplyr::summarise(n_loci=sum(Reads >= 100)) %>% pivot_wider(names_from = Pool, values_from = n_loci) %>% data.frame()
colnames(loci_stats)=c("SampleID",pool_columns)


df1=read.delim(samplestatFILE,header=T)
sample_stats=df1 %>% 
  pivot_wider(names_from = Stage, values_from = Reads) %>%
  data.frame() %>%
  mutate(SampleNumber=sapply(str_split(SampleID,'_S'),tail,1)) %>%
  mutate(SampleID=sapply(str_split(SampleID,'_S(\\d+)'),head,1))
sample_stats$SampleID=factor(sample_stats$SampleID,
                              levels=unique(samples$SampleID)) 
sample_stats=sample_stats[,c("SampleNumber","SampleID","Input","No.Dimers", "Amplicons")]%>% 
  arrange(SampleNumber)

colnames(sample_stats) = c("#","Sample","Input","No Dimers","Amplicons")

ampdata=read.delim(ampliconFILE,header=T)

#Histogram#
p1=ggplot(data=df, aes(x=Reads+0.1)) +  
  geom_histogram( ) + 
  scale_y_continuous() +
  scale_x_log10()+
  guides(fill=FALSE) + 
  xlab("Read Count") + 
  ylab("Frequency") + 
  ggtitle("\nNumber of Reads/Locus") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~SampleID,ncol=6) + 
  theme(plot.title = element_text(hjust = 0.5, size=25)) + 
  theme(strip.text.x = element_text(size = 25)) +
  theme(axis.text.x = element_text(size = 25)) + 
  theme(axis.text.y = element_text(size = 25)) + 
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.title.y = element_text(size = 25))+
  scale_y_log10()

ggsave(file="quality_report/reads_histograms.pdf", width=40, height=60, dpi=300, limitsize = FALSE)


#Boxplot#
numsamples=length(unique(df$SampleID))
samples_group1=unique(df$SampleID)[1:ceiling(numsamples/2)]
samples_group2=unique(df$SampleID)[(ceiling(numsamples/2)+1):numsamples]

p2a=ggplot( data=df %>% dplyr::filter(SampleID %in% samples_group1),aes(x=SampleID, y=Reads)) +
    geom_boxplot(color="#993333",lwd=0.75) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none",plot.title = element_text(size=11)) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Number of Reads/Locus") +
    xlab("") + facet_wrap(~Pool,ncol=1)
    
p2b=ggplot( data=df %>% dplyr::filter(SampleID %in% samples_group2),aes(x=SampleID, y=Reads)) +
    geom_boxplot(color="#993333",lwd=0.75) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none",plot.title = element_text(size=11)) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Number of Reads/Locus") +
    xlab("") + facet_wrap(~Pool,ncol=1)  
    
   
df2=df
df2$Reads[which(df$Reads == 0)]=0.1  
p3=ggplot(df2) +   
  ggbeeswarm::geom_quasirandom(aes(x=1,y=Reads,color = Pool),dodge.width = 0.5,size=3)+
  scale_y_log10()+
  facet_wrap(~SampleID,ncol=6)+
  theme_bw() +
  xlab("")+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
  geom_hline(yintercept = 100,linetype="dashed",color = "grey")+
  theme(strip.text.x = element_text(size = 25))+
  theme(axis.text.y = element_text(size = 25)) + 
  theme(axis.title.y = element_text(size = 25))+
  theme(plot.title = element_text(hjust = 0.5, size=30)) + 
  theme(legend.text = element_text(size = 25), 
        legend.title = element_text(size = 25), 
        legend.position="bottom")

ggsave(file="quality_report/swarm_plots.pdf", width=60, height=160, dpi=300, limitsize=FALSE)

#Length vs. Reads#
df1=df %>% left_join(ampdata,by = c("Locus" = "amplicon")) %>% select(SampleID,Locus,Reads,ampInsert_length,Pool) %>% data.frame()
p4=ggplot(df1,aes(x=ampInsert_length,y=Reads+0.1,color = Pool)) + ggtitle("Locus Length vs. Reads") + 
  geom_point(alpha=0.9,size=2.5) + 
  scale_y_log10()+
  xlab("Locus Insert Length") + 
  theme_bw() + 
  facet_wrap(~SampleID,ncol=6) + 
  theme(strip.text.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 25)) + 
  theme(axis.text.y = element_text(size = 25)) + 
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.title.y = element_text(size = 20))+
  theme(plot.title = element_text(hjust = 0.5, size=30)) + 
  theme(legend.text = element_text(size = 25), 
        legend.title = element_text(size = 25), 
        legend.position="bottom") 

ggsave(file="quality_report/length_vs_reads.pdf", width=60, height=200, dpi=300, limitsize=FALSE)

#pdf(paste(outDIR,"/QCplots.pdf",sep=""), onefile = TRUE)
#grid.arrange(p1,p2a,p2b,p3,p4)
#dev.off()

#p=list(p1,p3,p4)
#ml <- marrangeGrob(p, nrow=1, ncol=1)
## non-interactive use, multipage pdf
#ggsave(filename = paste(outDIR,"/QCplots.pdf",sep=""), plot = ml, width = 15, height = 9)




##RMarkdown report##

currentDate <- Sys.Date()
#rmd_file=paste(outDIR,"/QCplots.Rmd",sep="")
rmd_file=paste(outDIR,"QCplots.Rmd",sep='/')
file.create(rmd_file)
p=list(sample_stats,sample_amplicon_stats,loci_stats,p1,p3,p4)

file <- tempfile()
saveRDS(p, file)


c(paste0("---\ntitle: \"QC summary Report\"\nauthor: \"GenMoz\"\ndate: ",
currentDate,
"\noutput: html_document\n---\n")) %>% write_lines(rmd_file)
#c("```{r echo=FALSE, message=FALSE, warning=FALSE}\nplot_list=readRDS(file)\nlapply(plot_list,print)\n```") %>% write_lines(rmd_file,append=T)
c("```{r echo=FALSE, results=\'asis\', message=FALSE, warning=FALSE}\nplot_list=readRDS(file)\nknitr::kable(plot_list[1], caption=\"Cutadapt Sample Statistics\")\n```") %>% write_lines(rmd_file,append=T)
c("```{r echo=FALSE, message=FALSE, warning=FALSE}\nplot_list=readRDS(file)\nknitr::kable(plot_list[2], caption=\"Sample Median Reads per Pool\")\n```") %>% write_lines(rmd_file,append=T)
c("```{r echo=FALSE, message=FALSE, warning=FALSE}\nplot_list=readRDS(file)\nknitr::kable(plot_list[3], caption=\"Sample Number of Loci with 100 Reads or More per Pool\")\n```") %>% write_lines(rmd_file,append=T)
c("```{r echo=FALSE, message=FALSE, warning=FALSE}\nplot_list[-c(1:3)]\n```") %>% write_lines(rmd_file,append=T)

    rmarkdown::render(rmd_file, params = list(file=file, output_file = html_document()))
