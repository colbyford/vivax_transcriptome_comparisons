## Top gene selection using Schizont cell type percentages

library(readr)
library(dplyr)
# BiocManager::install("edgeR")
library(edgeR)

## Example: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

## Load in Data


### Ethiopian Isolates
cibersort_results <- read_csv('output/CIBERSORTx_Job3_Results.csv')
gene_counts <- read.delim('output/feature_counts.csv', header = TRUE, sep = ",", row.names = "gene_id")

### Cambodian Isolates
cibersort_results <- read_csv('output/CIBERSORTx_Job4_Results.csv')
gene_counts <- read.delim('output/feature_counts.csv', header = TRUE, sep = ",", row.names = "gene_id")

### Maryland Set Isolates (Tebben, et al.)
# cibersort_results <- read_csv('one-to-one_test/pb_orthologsonly.txt/CIBERSORTx_Job15_Results.csv')
# gene_counts <- read.delim('one-to-one_test/pb_orthologsonly.txt/pb_orthologsonly_pvgenes.txt', header = TRUE, sep = "\t", row.names = "Gene")

## Convert data into DGEList
y <- DGEList(counts=gene_counts)
y <- calcNormFactors(y)

##  Estimating dispersions
y <- estimateCommonDisp(y)

## 4.1.6 Design Matrix (use Schizont percentages from CIBERSORTx)
design <- model.matrix(~cibersort_results$Schizonts)
rownames(design) <- colnames(y)

## 4.1.7 Estimate Dispersion
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)


## 4.1.8 Differential Expression
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt, n=30)


fit_outputs <- lrt$table %>% 
  mutate(gene_id = row.names(lrt$table),
         coeff = lrt$coefficients[,2]) %>% 
  filter(PValue <= 0.05)

o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]


summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-1, 1), col="blue")
