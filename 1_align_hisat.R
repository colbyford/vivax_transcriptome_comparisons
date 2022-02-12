library(Rhisat2)

## Building a genome index
refs <- list.files("refs", full.names=TRUE, pattern="\\.fasta$")

dir <- "output"

hisat2_build(references=refs, outdir=dir, prefix="vivax_index", 
             force=TRUE, strict=TRUE, execute=TRUE)

## Aligning reads to the genome index
samples <- read.table("sample_summary.txt", header = TRUE, stringsAsFactors = FALSE)

for (i in seq_along(samples$SampleID)){
  cat("Processing Sample:", samples$SampleID[i], "\n")

  reads <- as.list(samples[i, 2:3])
  
  hisat2(sequences=reads,
         index=file.path(dir, "vivax_index"),
         type="paired",
         outfile=file.path(dir, paste0(samples$SampleID[i], "_output.sam")),
         threads = parallel::detectCores()-1,
         force=TRUE, strict=TRUE, execute=TRUE)
}
  