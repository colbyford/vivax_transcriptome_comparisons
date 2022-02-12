library(Rsamtools)
library(stringr)

## List all output SAM files
sam_files <- list.files("output", full.names=TRUE, pattern="\\.sam$")

## Convert all SAM files to BAM
cat("Converting", length(sam_files), "SAM files to BAM.\n")

lapply(sam_files, asBam)


## List all output BAM files
bam_files <- list.files("output", full.names=TRUE, pattern="\\.bam$")


## Filter down to only P. vivax reads (MAPQ > 30)
filter <- FilterRules(list(MinQuality = function(x) x$mapq > 30))
# filter <- FilterRules(list(MinQuality = function(x) x$mapq > 0))

## Initialize filtering statistics dataframe
filter_stats <- data.frame(file_original = character(0L),
                           file_filtered = character(0L),
                           records_in_original = integer(0L),
                           records_in_filtered = integer(0L))

for (i in seq_along(bam_files)){
  cat("Processing BAM file:", bam_files[i], "\n")
  
  bam_iter <- bam_files[i]
  
  ## Define output file location/name
  desitination_output <- paste0(str_remove(bam_iter, ".bam"), "_filtered.bam")
  
  
  ## Perform Filtering
  filterBam(bam_iter,
            destination = desitination_output,
            filter = filter)
  
  ## Get record counts before and after filtering
  count_original <- countBam(bam_iter)
  count_filtered <- countBam(desitination_output)
  
  ## Append to filtering statistics dataframe
  filter_stats_iter <- data.frame(file_original = count_original$file,
                                  file_filtered = count_filtered$file,
                                  records_in_original = count_original$records,
                                  records_in_filtered = count_filtered$records)
  
  filter_stats <- rbind(filter_stats, filter_stats_iter)
  
}

## Write out filtering statistics dataframe
write.csv(filter_stats, "output/filter_stats.csv", row.names = FALSE)
