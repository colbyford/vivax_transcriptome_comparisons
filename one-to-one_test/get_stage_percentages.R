## Derive Stage-level expression by gene

library(dplyr)
library(readxl)
library(tidyr)

pd_ortho <- read_excel("pb_orthologsonly.txt/pb_orthologsonly_transpose.xlsx")

gene_ids <- colnames(pd_ortho)

pd_ortho_pvt <- pd_ortho %>% pivot_longer(!Stage, names_to = "gene_id")

pd_ortho_pvt_wide <- pd_ortho_pvt %>%
  pivot_wider(names_from = "Stage",
              values_from = "value",
              values_fn=sum) %>% 
  select(!c("F", "M")) %>% 
  rowwise() %>% 
  mutate(total = sum(c(troph, ring, schiz)),
         troph_pct = troph/total,
         ring_pct = ring/total,
         schiz_pct = schiz/total)

readr::write_csv(pd_ortho_pvt_wide, "stage_percentages_by_gene.csv")



## t-SNE
# BiocManager::install("M3C")
library(M3C)
data(mydata)

pd_ortho_orig <- read.delim("pb_orthologsonly.txt/pb_orthologsonly.txt", sep = "\t")

stages <- pd_ortho$Stage

tsne(pd_ortho_orig %>% select(!Gene),
     labels = stages)


