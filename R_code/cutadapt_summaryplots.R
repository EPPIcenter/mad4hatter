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



df=read.table(summaryFILE,header=T)
df$SampleName=as.factor(df$SampleName)
df$Amplicon=as.factor(df$Amplicon)
df = df %>% 
  mutate(SampleNumber=sapply(str_split(SampleName,'_S'),tail,1)) %>% 
  mutate(SampleName=sapply(str_split(SampleName,'_S(\\d+)'),head,1)) %>% 
  mutate(Pool=sapply(str_split(Amplicon,'-'),tail,1)) %>% 
  arrange(SampleName) %>% 
  data.frame()

samples = df %>% select(SampleName,SampleNumber) %>% distinct() %>% 
  arrange(SampleNumber)

df$SampleName=factor(df$SampleName,levels=unique(samples$SampleName))

amplicon_stats=df %>% select(-Pool) %>% pivot_wider(names_from = Amplicon, values_from = NumReads) %>% data.frame()
write.table(amplicon_stats, file=paste(outDIR,"/amplicon_stats.txt",sep=""), quote=F, sep ="\t", col.names=T, row.names=F)

sample_amplicon_stats=df %>% group_by(SampleName,Pool) %>% dplyr::summarise(medianReads=median(NumReads)) %>% pivot_wider(names_from = Pool, values_from = medianReads) %>% data.frame()

loci_stats= df %>% group_by(SampleName) %>% dplyr::summarise(percent_loci_with_100_reads = (sum(NumReads >= 100) / n()) * 100) %>% data.frame()

colnames(sample_amplicon_stats)=c("SampleName","Pool_1A","Pool_1AB","Pool_1B","Pool_1B2", "Pool_2")

df1=read.delim(samplestatFILE,header=T)
sample_stats=df1 %>% 
  pivot_wider(names_from = X, values_from = NumReads) %>%
  data.frame() %>%
  mutate(SampleNumber=sapply(str_split(SampleName,'_S'),tail,1)) %>%
  mutate(SampleName=sapply(str_split(SampleName,'_S(\\d+)'),head,1))
sample_stats$SampleName=factor(sample_stats$SampleName,
                              levels=unique(samples$SampleName)) 
sample_stats=sample_stats[,c("SampleNumber","SampleName","Input","No.Dimers","Amplicons")]%>% 
  arrange(SampleNumber)

colnames(sample_stats) = c("#","Sample","Input","No Dimers","Amplicons")

ampdata=read.delim(ampliconFILE,header=T)

#Histogram#
p1=ggplot(data=df, aes(x=NumReads+0.1)) +  
  geom_histogram( ) + 
  scale_y_continuous() +
  scale_x_log10()+
  guides(fill=FALSE) + 
  xlab("Read Count") + 
  ylab("Frequency") + 
  ggtitle("\nNumber of Reads/Amplicon") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~SampleName,ncol=6) + 
  theme(plot.title = element_text(hjust = 0.5, size=25)) + 
  theme(strip.text.x = element_text(size = 25)) +
  theme(axis.text.x = element_text(size = 25)) + 
  theme(axis.text.y = element_text(size = 25)) + 
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.title.y = element_text(size = 25))+
  scale_y_log10()

ggsave(file="quality_report/reads_histograms.pdf", width=40, height=60, dpi=300, limitsize = FALSE)


#Boxplot#
numsamples=length(unique(df$SampleName))
samples_group1=unique(df$SampleName)[1:ceiling(numsamples/2)]
samples_group2=unique(df$SampleName)[(ceiling(numsamples/2)+1):numsamples]

p2a=ggplot( data=df %>% dplyr::filter(SampleName %in% samples_group1),aes(x=SampleName, y=NumReads)) +
    geom_boxplot(color="#993333",lwd=0.75) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none",plot.title = element_text(size=11)) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Number of Reads/Amplicon") +
    xlab("") + facet_wrap(~Pool,ncol=1)
    
p2b=ggplot( data=df %>% dplyr::filter(SampleName %in% samples_group2),aes(x=SampleName, y=NumReads)) +
    geom_boxplot(color="#993333",lwd=0.75) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none",plot.title = element_text(size=11)) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Number of Reads/Amplicon") +
    xlab("") + facet_wrap(~Pool,ncol=1)  
    
   
df2=df
df2$NumReads[which(df$NumReads == 0)]=0.1  
p3=ggplot(df2) +   
  ggbeeswarm::geom_quasirandom(aes(x=1,y=NumReads,color = Pool),dodge.width = 0.5,size=3)+
  scale_y_log10()+
  facet_wrap(~SampleName,ncol=6)+
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

#Length vs. NumReads#
df1=df %>% left_join(ampdata,by = c("Amplicon" = "amplicon")) %>% select(SampleName,Amplicon,NumReads,ampInsert_length,Pool) %>% data.frame()
p4=ggplot(df1,aes(x=ampInsert_length,y=NumReads+0.1,color = Pool)) + ggtitle("Amplicon Length vs. NumReads") + 
  geom_point(alpha=0.9,size=2.5) + 
  scale_y_log10()+
  xlab("Amplicon Insert Length") + 
  theme_bw() + 
  facet_wrap(~SampleName,ncol=6) + 
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
c("```{r echo=FALSE, message=FALSE, warning=FALSE}\nplot_list=readRDS(file)\nknitr::kable(plot_list[3], caption=\"Sample Percent of Loci with 100 Reads or More\")\n```") %>% write_lines(rmd_file,append=T)
c("```{r echo=FALSE, message=FALSE, warning=FALSE}\nplot_list[-c(1:3)]\n```") %>% write_lines(rmd_file,append=T)

    rmarkdown::render(rmd_file, params = list(file=file, output_file = html_document()))
