demuxed_folder = "~/Data/test_times/results/cutadapt_TES_quarter_demux"

fastq.summary = ShortRead::countFastq(demuxed_folder,pattern = "R1")

fastq.reads = data.frame(filename = rownames(fastq.summary),
                            reads = fastq.summary$records)

fastq.reads$locus = sapply(strsplit(fastq.reads$filename,"_trimmed"),"[",1)

fastq.reads = fastq.reads[,c("locus","reads")]

write.table(fastq.reads, file='~/Data/test_times/results/cutadapt_TES_quarter_demux/count_reads.tsv', quote=FALSE, sep='\t')

