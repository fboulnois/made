#' Microarray pipeline
#'
#' Run the full microarray pipeline.
#'
#' Runs the entire microarray pipeline according to the options set in the
#' config. The pipeline consists of quality assessment, normalization,
#' summarization, and creation of the final report.
#'
#' @param config Character string consisting of the path to the configuration
#' file generated using the \code{write.yaml.config} function or parsed
#' configuration list associated with a microarray experiment.
#'
#' @seealso \code{\link{write.yaml.config}} to generate the configuration file,
#' \code{\link{ma.normalize}} to perform normalization,
#' \code{\link{ma.summarize}} to perform gene level summarization, and
#' \code{\link{write.report}} to generate the final report including the quality
#' assessment.
#'
#' @export
ma.pipeline <- function(config)
{
  config <- read.yaml.config(config)

  eset <- ma.normalize(config)

  if(config$pipeline$summarization)
  {
    results <- ma.summarize(config, eset)
  }

  if(config$pipeline$write_report)
  {
    report <- write.report(config, eset, results)
  }

  return(report)
}
