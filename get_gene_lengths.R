## Transcript Lengths

library(GenomicFeatures)
library(readr)


gff <- "refs/PlasmoDB-48_PvivaxP01.gff"

txdb <- makeTxDbFromGFF(
  file = gff,
  dataSource = "PlasmoDB 48",
  organism = "Plasmodium vivax"
)

tx_by_gene <- transcriptsBy(txdb, by="gene")
gene_lens <- max(width(tx_by_gene))

genes <- genes(txdb)
gene_lengths <- width(genes)


genes_df <- as.data.frame(genes)

write_csv(genes_df, "refs/gene_lengths.csv")

