## Set CRAN Mirror
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)


## Install Packages
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

install.packages("stringr")
install.packages("pheatmap")


# ShortRead: http://bioconductor.org/packages/release/bioc/vignettes/ShortRead/inst/doc/Overview.pdf
BiocManager::install("ShortRead")

# Hisat2: https://bioconductor.org/packages/release/bioc/vignettes/Rhisat2/inst/doc/Rhisat2.html
# devtools::install_github("fmicompbio/Rhisat2")
BiocManager::install("Rhisat2")

# Rsamtools: https://bioconductor.org/packages/release/bioc/vignettes/Rsamtools/inst/doc/Rsamtools-Overview.pdf
BiocManager::install("Rsamtools")

# Rsubread: https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf
BiocManager::install("Rsubread")

# DESeq2: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# devtools::install_github("mikelove/DESeq2")
BiocManager::install("DESeq2")



