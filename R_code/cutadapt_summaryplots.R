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


df=read.table(summaryFILE,header=F)
colnames(df)=c("SampleName","Amplicon","NumReads")
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

df1=read.delim(samplestatFILE,header=F)
sample_stats=df1 %>% 
  pivot_wider(names_from = V2, values_from = V3) %>% 
  dplyr::rename(SampleName = V1) %>% data.frame() %>% 
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
  ggtitle("\nNumber of reads/Amplicon") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 6)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~SampleNumber,ncol=6) + 
  theme(strip.text.x = element_text(size = 7))


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
  ggbeeswarm::geom_quasirandom(aes(x=1,y=NumReads,color = Pool),dodge.width = 0.5,size=0.5)+
  scale_y_log10()+
  facet_wrap(~SampleNumber,ncol=8)+
  theme_bw() +
  xlab("")+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
  geom_hline(yintercept = 100,linetype="dashed",color = "grey")+
  theme(strip.text.x = element_text(size = 7))+
  theme(legend.position = "bottom")


#Length vs. NumReads#
df1=df %>% left_join(ampdata,by = c("Amplicon" = "amplicon")) %>% select(SampleNumber,Amplicon,NumReads,ampInsert_length,Pool) %>% data.frame()
p4=ggplot(df1,aes(x=ampInsert_length,y=NumReads+0.1,color = Pool)) + ggtitle("Amplicon Length vs. NumReads") + 
  geom_point(alpha=0.9,size=0.5) + 
  scale_y_log10()+
  xlab("Amplicon Insert Length") + 
  theme_bw() + 
  theme(axis.text.y  = element_text(size = 7)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~SampleNumber,ncol=8) + 
  theme(strip.text.x = element_text(size = 7))+
  theme(legend.position = "bottom")

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
p=list(sample_stats,p1,p3,p4)

file <- tempfile()
saveRDS(p, file)


c(paste0("---\ntitle: \"QC summary Report\"\nauthor: \"GenMoz\"\ndate: ",
currentDate,
"\noutput: html_document\n---\n")) %>% write_lines(rmd_file)
#c("```{r echo=FALSE, message=FALSE, warning=FALSE}\nplot_list=readRDS(file)\nlapply(plot_list,print)\n```") %>% write_lines(rmd_file,append=T)
c("```{r echo=FALSE, results=\'asis\', message=FALSE, warning=FALSE}\nplot_list=readRDS(file)\nknitr::kable(plot_list[1], caption=\"Cutadapt Sample Statistics\")\n```") %>% write_lines(rmd_file,append=T)
c("```{r echo=FALSE, message=FALSE, warning=FALSE}\nplot_list[c(2:4)]\n```") %>% write_lines(rmd_file,append=T)

    rmarkdown::render(rmd_file, params = list(file=file, output_file = html_document()))