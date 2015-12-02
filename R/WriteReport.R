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
#' @param results A summary of all the data in the expression set for each
#' group comparison generated using the \code{ma.summarize} function.
#'
#' @importFrom rmarkdown render
#'
#' @export
write.report <- function(config, eset, results)
{
  config <- read.yaml.config(config)

  if(missing(eset) || class(eset) != "ExpressionSet")
  {
    stop("Expression set is required to write report.")
  }

  if(missing(results) || class(results) != "MadeSummary")
  {
    stop("Experiment summary is required to write report.")
  }

  infile  <- normalizePath(system.file("extdata", "report-doc.rmd", package = "made"), mustWork = TRUE)
  outfile <- normalizePath(file.path(dirname(config$groups$group_file), sprintf("%s_microarray_report.html", .timestamped())), mustWork = FALSE)

  tryCatch({
    rmarkdown::render(input = infile, output_file = outfile, clean = FALSE)
  }, error = function(e)
  {
    emsg <- e$message
    if(grepl("pandoc .* not found", emsg))
    {
      emsg <- paste(emsg, "Download it from https://github.com/jgm/pandoc/releases")
    }
    stop(emsg)
  })

  return(outfile)
}
