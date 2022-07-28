library(tidyverse)
library(gridExtra)
library(ggbeeswarm)
library(rmarkdown)
library(knitr)
library(ggforce)

args = commandArgs(trailingOnly=T)

summaryFILE=args[1]
samplestatFILE=args[2]
ampliconFILE=args[3]
outDIR=args[4]


df=read.delim(summaryFILE,header=F)
colnames(df)=c("SampleName","Amplicon","NumReads")
df$SampleName=as.factor(df$SampleName)
df$Amplicon=as.factor(df$Amplicon)
df = df %>% mutate(SampleName=sapply(str_split(SampleName,'_'),head,1)) %>% mutate(Pool=sapply(str_split(Amplicon,'-'),tail,1)) %>% arrange(SampleName) %>% data.frame()
df$SampleName=factor(df$SampleName,levels=unique(df$SampleName))

amplicon_stats=df %>% select(-Pool) %>% pivot_wider(names_from = Amplicon, values_from = NumReads) %>% data.frame()
write.table(amplicon_stats, file=paste(outDIR,"/amplicon_stats.txt",sep=""), quote=F, sep ="\t", col.names=T, row.names=F)

sample_amplicon_stats=df %>% group_by(SampleName,Pool) %>% dplyr::summarise(medianReads=median(NumReads)) %>% pivot_wider(names_from = Pool, values_from = medianReads) %>% data.frame()
colnames(sample_amplicon_stats)=c("SampleName","Pool_1A","Pool_1AB","Pool_1B","Pool_2")

df1=read.delim(samplestatFILE,header=F)
sample_stats=df1 %>% pivot_wider(names_from = V2, values_from = V3) %>% dplyr::rename(SampleName = V1) %>% data.frame()


ampdata=read.delim(ampliconFILE,header=T)

#setwd(outDIR)
numsamples=length(unique(df$SampleName))
numcols=8
numrows=ifelse(numsamples>=48,6,ceiling(numsamples/numcols))

numpages=ifelse(numsamples>=144,4,ifelse(numsamples>=96,3,ifelse(numsamples>=48,2,1)))

#numcols=ifelse(numsamples>=32,12,ifelse(numsamples>=16,8,4))


#Histogram#
#p1=ggplot(data=df, aes(x=NumReads)) +  geom_histogram( binwidth=100, fill="#993333", color="#990000", alpha=0.9) + scale_y_continuous(breaks=seq(0,10,5),limits=c(0,10)) +  guides(fill=FALSE) + xlab("Read Count") + ylab("Frequency") + ggtitle("\nNumber of reads/Amplicon") + theme_bw() + theme(axis.text.x = element_text(size = 6)) + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~SampleName,ncol=numcols) + theme(strip.text.x = element_text(size = 7)) +   theme(strip.background = element_blank(),strip.text = element_text(size = rel(0.8), margin = margin()),panel.spacing = unit(3, "pt"))

p1=ggplot(data=df, aes(x=NumReads)) +  geom_histogram( binwidth=100, fill="#993333", color="#990000", alpha=0.9) + scale_y_continuous(breaks=seq(0,10,5),limits=c(0,10)) +  guides(fill=FALSE) + xlab("Read Count") + ylab("Frequency") + ggtitle("\nNumber of reads/Amplicon") + theme_bw() + theme(axis.text.x = element_text(size = 6)) + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap_paginate(~SampleName,ncol=numcols,nrow=numrows) + theme(strip.text.x = element_text(size = 7)) +   theme(strip.background = element_blank(),strip.text = element_text(size = rel(0.8), margin = margin()),
panel.spacing = unit(3, "pt"))

p1L <- vector(mode = "list", length = numpages)

for(kk in 1:numpages)
{
p1L[[kk]]=p1+facet_wrap_paginate(~SampleName,ncol=numcols,nrow=numrows,page=kk)
}

#p11=p1+facet_wrap_paginate(~SampleName,ncol=8,nrow=6,page=1)
#p12=p1+facet_wrap_paginate(~SampleName,ncol=8,nrow=6,page=2)
#p13=p1+facet_wrap_paginate(~SampleName,ncol=8,nrow=6,page=3)
#p14=p1+facet_wrap_paginate(~SampleName,ncol=8,nrow=6,page=4) 

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
    
   
 df$NumReads[which(df$NumReads == 0)]=0.1  
 #p3=ggplot(df) + ggbeeswarm::geom_quasirandom(aes(x=1,y=NumReads,color = Pool),dodge.width = 0.5,size=0.5)+ scale_y_log10()+ facet_wrap(~SampleName,ncol=numcols)+ xlab("")+ theme_bw() + theme(strip.background = element_blank(),strip.text = element_text(size = rel(0.7), margin = margin()),panel.spacing = unit(3, "pt")) + geom_hline(yintercept = 100,linetype="dashed",color = "grey")
      
p3=ggplot(df) + ggbeeswarm::geom_quasirandom(aes(x=1,y=NumReads,color = Pool),dodge.width = 0.5,size=0.5)+ scale_y_log10()+ facet_wrap_paginate(~SampleName,ncol=numcols,nrow=numrows) + xlab("")+ theme_bw() + theme(strip.background = element_blank(),strip.text = element_text(size = rel(0.7), margin = margin()),panel.spacing = unit(3, "pt")) + geom_hline(yintercept = 100,linetype="dashed",color = "grey")

p3L <- vector(mode = "list", length = numpages)

for(kk in 1:numpages)
{
p3L[[kk]]=p3+facet_wrap_paginate(~SampleName,ncol=numcols,nrow=numrows,page=kk)
}

      
#p31=p3+facet_wrap_paginate(~SampleName,ncol=8,nrow=6,page=1)
#p32=p3+facet_wrap_paginate(~SampleName,ncol=8,nrow=6,page=2)
#p33=p3+facet_wrap_paginate(~SampleName,ncol=8,nrow=6,page=3)
#p34=p3+facet_wrap_paginate(~SampleName,ncol=8,nrow=6,page=4) 

#Length vs. NumReads#
df1=df %>% left_join(ampdata,by = c("Amplicon" = "amplicon")) %>% select(SampleName,Amplicon,NumReads,ampInsert_length) %>% data.frame()

#p4=ggplot(df1,aes(x=ampInsert_length,y=NumReads)) + ggtitle("Amplicon Length vs. NumReads") + geom_point(alpha=0.9,color="#993333") + xlab("Amplicon Insert Length") + theme_bw() + theme(axis.text.y  = element_text(size = 7)) + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~SampleName,ncol=numcols) + theme(strip.text.x = element_text(size = 7)) + theme(strip.background = element_blank(),strip.text = element_text(size = rel(0.8), margin = margin()),panel.spacing = unit(3, "pt"))

p4=ggplot(df1,aes(x=ampInsert_length,y=NumReads)) + ggtitle("Amplicon Length vs. NumReads") + geom_point(alpha=0.9,color="#993333") + xlab("Amplicon Insert Length") + theme_bw() + theme(axis.text.y  = element_text(size = 7)) + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap_paginate(~SampleName,ncol=numcols,nrow=numrows) + theme(strip.text.x = element_text(size = 7)) + theme(strip.background = element_blank(),strip.text = element_text(size = rel(0.8), margin = margin()),panel.spacing = unit(3, "pt"))


p4L <- vector(mode = "list", length = numpages)

for(kk in 1:numpages)
{
p4L[[kk]]=p4+facet_wrap_paginate(~SampleName,ncol=numcols,nrow=numrows,page=kk)
}

#p41=p4+facet_wrap_paginate(~SampleName,ncol=8,nrow=6,page=1)
#p42=p4+facet_wrap_paginate(~SampleName,ncol=8,nrow=6,page=2)
#p43=p4+facet_wrap_paginate(~SampleName,ncol=8,nrow=6,page=3)
#p44=p4+facet_wrap_paginate(~SampleName,ncol=8,nrow=6,page=4) 

#pdf(paste(outDIR,"/QCplots.pdf",sep=""), onefile = TRUE)
#grid.arrange(p1,p2a,p2b,p3,p4)
#dev.off()

#p=list(p1,p3,p4)
#ml <- marrangeGrob(p, nrow=1, ncol=1)
## non-interactive use, multipage pdf
#ggsave(filename = paste(outDIR,"/QCplots.pdf",sep=""), plot = ml, width = 15, height = 9)


##RMarkdown report##

currentDate <- Sys.Date()
rmd_file=paste(outDIR,"/QCplots.Rmd",sep="")
#rmd_file="QCplots.Rmd"
file.create(rmd_file)
#p=list(sample_stats,p1,p3,p4)
#p=list(sample_stats,p11,p12,p13,p14,p31,p32,p33,p34,p41,p42,p43,p44)
p=list(sample_stats,sample_amplicon_stats,p1L,p3L,p4L)

file <- tempfile()
saveRDS(p, file)


c(paste0("---\ntitle: \"QC summary Report\"\nauthor: \"GenMoz\"\ndate: ",currentDate,"\noutput: html_document\n---\n")) %>% write_lines(rmd_file)
#c("```{r echo=FALSE, message=FALSE, warning=FALSE}\nplot_list=readRDS(file)\nlapply(plot_list,print)\n```") %>% write_lines(rmd_file,append=T)
c("```{r echo=FALSE, results=\'asis\', message=FALSE, warning=FALSE}\nplot_list=readRDS(file)\nknitr::kable(plot_list[1], caption=\"Cutadapt Sample Reads\")\n```") %>% write_lines(rmd_file,append=T)
c("```{r echo=FALSE, results=\'asis\', message=FALSE, warning=FALSE}\nplot_list=readRDS(file)\nknitr::kable(plot_list[2], caption=\"Sample Median Reads per Pool\")\n```") %>% write_lines(rmd_file,append=T)
c("```{r echo=FALSE, message=FALSE, warning=FALSE}\nplot_list[-c(1,2)]\n```") %>% write_lines(rmd_file,append=T)

    rmarkdown::render(rmd_file, params = list(file=file, output_file = html_document()))

