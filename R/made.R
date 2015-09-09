#' made: A package for the reproducible analysis of microarray data.
#'
#' The \code{\strong{made}} package provides tools for the reproducible analysis
#' of oligonucleotide arrays using a simple, clear, and consistent interface.
#'
#' A user should first specify the options they require and the groups or
#' conditions they want to compare and generate a configuration file using the
#' \code{\link{write.yaml.config}} function. After this they can run the entire
#' pipeline using the \code{\link{ma.pipeline}} function or integrate portions
#' of the analysis separately in a larger workflow. The pipeline includes
#' functions such as normalization using \code{\link{ma.normalize}}, gene level
#' summarization using \code{\link{ma.summarize}}, and generating a high-quality
#' report using \code{\link{write.report}}, including quality assessment if
#' desired. Other functions can be used to read and write various files required
#' by the analysis programatically. Reading the configuration or group file can
#' be accomplished using the functions \code{\link{read.yaml.config}} and
#' \code{\link{read.group.file}} respectively, while writing the configuration
#' file can be accomplished using the function \code{\link{write.yaml.config}}.
#'
#' @docType package
#' @name made
NULL
