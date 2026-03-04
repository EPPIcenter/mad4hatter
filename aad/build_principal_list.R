library(tidyverse)

principal_list = read.delim("~/Documents/repos/mad4hatter/panel_information/principal_resistance_marker_info_table.tsv")

resmarkers_aad = read.delim("~/Downloads/resistance_markers_amplicon_v4.txt") %>% 
  mutate(GeneID=paste0(strrep("0",7-nchar(as.character(GeneID))),
                      as.character(GeneID))
         )%>% 
  select(gene_id = GeneID,
         gene = Gene,
         aa_position = CodonID, 
         chrom = chr,
         start,
         stop, 
         strand) %>% 
  distinct()

new_principal_list = rbind(
  principal_list,
  resmarkers_aad %>% 
    left_join(df_resmarker %>% 
                mutate(a = 1),
              by = join_by(gene_id, gene, aa_position, chrom, start, stop, strand)
              ) %>% 
    filter(is.na(a)) %>% 
    select(-a)
)


write_tsv(new_principal_list,"~/Documents/repos/mad4hatter/panel_information/principal_resistance_marker_info_table.tsv")
