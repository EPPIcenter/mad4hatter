library(tidyverse)

files <- list.files(
  path = "Documents/repos/mad4hatter/",
  pattern = "^target_id_conversion_table\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)%>%
  # exclude any path containing /aad/
  .[!grepl("/aad/", .)]


combined <- map_dfr(files, function(f) {
  
  # Extract subfolder name (immediate parent folder)
  pool_name <- basename(dirname(f))
  
  read_tsv(f) %>%
    mutate(pool = pool_name)
}) %>% 
  group_by(new_name,old_name) %>% 
  summarize(pool = paste(pool,collapse = ",")) %>% 
  write_tsv("~/Documents/repos/mad4hatter/aad/target_id_conversion_table.tsv")

