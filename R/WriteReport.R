#' Write microarray report
#'
#' Generate the final microarray report.
#'
#' Generates the microarray report summarizing the information in the
#' experiment. The information in the report depends on the config and the
#' packages available.
#'
#' @param config Character string consisting of the path to the configuration
#' file generated using the \code{write.yaml.config} function or parsed
#' configuration list associated with a microarray experiment.
#'
#' @param eset Expression set object describing microarray experiment at the
#' level of probes.
#'
#' @param tf A list of top tables for each comparison of interest and their
#' genes, associated log-fold changes, and statistical values.
write.report <- function(config, eset, tf)
{
  config <- read.yaml.config(config)

  if(missing(eset))
  {
    stop("Expression set is required to write report.")
  }

  if(missing(tf))
  {
    stop("List of gene top tables is required to write report.")
  }

  infile  <- normalizePath(system.file("extdata", "report.rmd", package = "made"), mustWork = TRUE)
  outfile <- normalizePath(file.path(dirname(config$groups$group_file), sprintf("%s_microarray_report.html", .timestamped())), mustWork = FALSE)
  rmarkdown::render(input = infile, output_file = outfile)
}
