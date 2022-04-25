library(readr)
library(dplyr)
library(tidyr)


pb_data <- read_tsv("pb_orthologsonly_keppleset_042022_TRANSPOSE.txt")

one_to_one <- read_tsv("../OnetoOnetoOne_Orthologs_PbPfPv.txt", col_names = c("PB", "PF", "PV"))


pb_data_pvt = pb_data %>% pivot_longer(starts_with("PF")) %>% 
  group_by(Stage, name) %>% 
  summarize(count = sum(value)) %>% 
  pivot_wider(id_cols = name, names_from = Stage, values_from = count) %>% 
  inner_join(one_to_one, by = c("name" = "PF")) %>% 
  rename(gene_id = PV, Schizonts = schiz, Trophozoites = troph, Male_Gametocytes = M, Female_Gametocytes = "F", Rings  = ring) %>% 
  mutate(included = 1) %>% 
  select(gene_id, included, Rings, Trophozoites, Schizonts)


write_tsv(pb_data_pvt, "pb_orthologsonly_keppleset_042022_TEMPLATE.txt")
