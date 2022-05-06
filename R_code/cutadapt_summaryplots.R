library(tidyverse)
library(gridExtra)
library(ggbeeswarm)

args = commandArgs(trailingOnly=T)

summaryFILE=args[1]
ampliconFILE=args[2]
outDIR=args[3]

df=read.delim(summaryFILE,header=F)
colnames(df)=c("SampleName","Amplicon","NumReads")
df$SampleName=as.factor(df$SampleName)
df$Amplicon=as.factor(df$Amplicon)
df = df %>% mutate(SampleName=sapply(str_split(SampleName,'_'),head,1)) %>% mutate(Pool=sapply(str_split(Amplicon,'-'),tail,1)) %>% arrange(SampleName) %>% data.frame()
df$SampleName=factor(df$SampleName,levels=unique(df$SampleName))

ampdata=read.delim(ampliconFILE,header=T)

#Histogram#
p1=ggplot(data=df, aes(x=NumReads)) +  geom_histogram( binwidth=100, fill="#993333", color="#990000", alpha=0.9) + scale_y_continuous(breaks=seq(0,10,5),limits=c(0,10)) +  guides(fill=FALSE) + xlab("Read Count") + ylab("Frequency") + ggtitle("\nNumber of reads/Amplicon") + theme_bw() + theme(axis.text.x = element_text(size = 6)) + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~SampleName,ncol=8) + theme(strip.text.x = element_text(size = 7))



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
 p3=ggplot(df) + ggbeeswarm::geom_quasirandom(aes(x=1,y=NumReads,color = Pool),dodge.width = 0.5,size=0.5)+
  scale_y_log10()+
  facet_wrap(~SampleName,ncol=12)+
  xlab("")+ theme_bw() +
  theme(strip.text.x = element_blank()) +
  geom_hline(yintercept = 100,linetype="dashed",color = "grey")
      

#Length vs. NumReads#
df1=df %>% left_join(ampdata,by = c("Amplicon" = "amplicon")) %>% select(SampleName,Amplicon,NumReads,ampInsert_length) %>% data.frame()
p4=ggplot(df1,aes(x=ampInsert_length,y=NumReads)) + ggtitle("Amplicon Length vs. NumReads") + 
  geom_point(alpha=0.9,color="#993333") + xlab("Amplicon Insert Length") + theme_bw() + theme(axis.text.y  = element_text(size = 7)) + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~SampleName,ncol=8) + theme(strip.text.x = element_text(size = 7))

#pdf(paste(outDIR,"/QCplots.pdf",sep=""), onefile = TRUE)
#grid.arrange(p1,p2a,p2b,p3,p4)
#dev.off()

p=list(p1,p3,p4)
ml <- marrangeGrob(p, nrow=1, ncol=1)
## non-interactive use, multipage pdf
ggsave(filename = paste(outDIR,"/QCplots.pdf",sep=""), plot = ml, width = 15, height = 9)


