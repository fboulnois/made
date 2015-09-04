write.report <- function(config, eset, tf)
{
  config <- read.yaml.config(config)

  if(missing(eset))
  {
    stop("Expression set is required to write report.")
  }

  if(missing(tf))
  {
    stop("List of top tables is required to write report.")
  }

  infile  <- normalizePath(system.file("extdata", "report.rmd", package = "made"), mustWork = TRUE)
  outfile <- normalizePath(file.path(dirname(config$groups$group_file), "microarray_report.html"), mustWork = FALSE)
  rmarkdown::render(input = infile, output_file = outfile)
}
