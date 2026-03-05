read_callable_file <- function( fname ) {
  mydf <- read.table( fname, header = T )
  sample <- gsub(".+/(.+).callableloci.autosomal.txt","\\1",fname)
  mydf$Sample <- sample
  return( mydf )
}

read_metrics_file <- function( fname ) {
  if ( grepl("wgs_metrics.txt", fname)) {
    sample <- gsub(".+/(.+).wgs_metrics.txt","\\1",fname)
  }
  if ( grepl("coverage_metrics",fname)) {
    sample <- gsub(".+/(.+).bwamem2.samtools.+","\\1",fname)
  }
  myfile <- readLines(fname)
  start_index <- grep("## METRICS CLASS", myfile)
  end_index <- grep("## HISTOGRAM", myfile)
  metrics_data <- myfile[(start_index + 1):(end_index - 1)]
  mydf <- read.table(text = paste(metrics_data, collapse = "\n"), header = TRUE, sep = "\t")
  mydf$Sample <- sample
  return( mydf )
}
