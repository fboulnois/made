write.report <- function(config, eset, tf)
{
  config <- read.yaml.config(config)

  infile  <- normalizePath(system.file("extdata", "report.rmd", package = "made"), mustWork = TRUE)
  outfile <- normalizePath(file.path(dirname(config$groups$group_file), "microarray_report.html"), mustWork = FALSE)
  rmarkdown::render(input = infile, output_file = outfile)
}
