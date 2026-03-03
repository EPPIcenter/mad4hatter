# ONLY RUN ON VERSION WITH NO POOL TAGS ON THE FILE NAMES OTHERWISE REWRITE

library(tidyverse)


files <- list.files(
  path = "Documents/repos/mad4hatter/",
  pattern = "^target_id_conversion_table\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)%>%
  # exclude any path containing /aad/
  .[!grepl("/aad/", .)]

# 2. Read and bind them, adding pool column
combined <- map_dfr(files, function(f) {
  
  # Extract subfolder name (immediate parent folder)
  pool_name <- basename(dirname(f))
  
  read_tsv(f) %>%
    mutate(pool = pool_name,
           row_id = row_number())
}) %>% 
  group_by(new_name,old_name) %>% 
  summarize(pool = paste(pool,collapse = ","),
            row_id = paste(row_id,collapse = ","))%>%
  mutate(
    old_name_fixed = case_when(
      
      # 1A cases 
      pool %in% c("D1.1",
                  "D1.1,M1.1",
                  "D1.1,M2.1") ~ str_c(old_name, "-1A"),
      
      # 1B cases 
      pool %in% c("R1.1",
                  "R1.1,R1.2",
                  "M1.1,R1.1,R1.2") ~ str_c(old_name, "-1B"),
      
      # 1AB cases 
      pool %in% c("D1.1,R1.1",
                  "D1.1,M1.1,R1.1",
                  "D1.1,M1.1,R1.1,R1.2") ~ str_c(old_name, "-1AB"),
      
      # 1B2 case 
      pool == "M1.1,R1.1,R1.2,R2.1" ~ str_c(old_name, "-1B2"),
      
      # -2 cases 
      pool %in% c("R2.1",
                  "M1.1,R2.1",
                  "M2.1,R2.1") ~ str_c(old_name, "-2"),
      
      # -6 cases 
      pool %in% c("M1.1",
                  "M1.addon",
                  "M2.1") ~ str_c(old_name, "-6"),
      TRUE ~ old_name
    )
  )




parent_dir <- "Documents/repos/mad4hatter/panel_information/"

# Expand comma-separated pools
expanded <- combined %>%
  separate_rows(pool,row_id, sep = ",")%>%
  mutate(pool = str_trim(pool),
         row_id = as.numeric(str_trim(row_id))
  )


# Get unique pool names
unique_pools <- unique(expanded$pool)

# Loop and write files
for (p in unique_pools) {
  
  df_pool <- expanded %>%
    filter(pool == p) %>% 
    arrange(row_id) %>%
    select(new_name, old_name=old_name_fixed)

  output_path <- file.path(
    parent_dir,
    p,
    "target_id_conversion_table.tsv"
  )
  
  write_tsv(df_pool, output_path)
}



