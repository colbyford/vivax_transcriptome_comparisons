library(Rsubread)

## List all output filtered BAM files
bam_files <- list.files("output", full.names=TRUE, pattern="\\_filtered.bam$")
# bam_files <- list.files("output", full.names=TRUE, pattern="\\_output.bam$")

## Count gene features from P01 annotation
feature_counts_iter <- featureCounts(bam_files,
                                     annot.ext="refs/PlasmoDB-48_PvivaxP01.gff",
                                     GTF.attrType = "ID",
                                     GTF.featureType = "gene",
                                     chrAliases = "refs/chr_aliases.csv",
                                     isPairedEnd = TRUE,
                                     isGTFAnnotationFile = TRUE,
                                     verbose = TRUE)

feature_counts_df <- feature_counts_iter$counts %>%
  # t() %>% 
  as.data.frame()

## Make sure gene_id is included (instead of as row names)
feature_counts_df$gene_id <- row.names(feature_counts_df)
row.names(feature_counts_df) <- NULL

## Move gene_id to the first column
feature_counts_df <- feature_counts_df[, c(ncol(feature_counts_df), 1:(ncol(feature_counts_df)-1))]


## Write out feature courts dataframe
write.csv(feature_counts_df, "output/feature_counts.csv", row.names = FALSE)
