#' Microarray normalization
#'
#' Background correct and normalize a microarray experiment.
#'
#' Reads in all the samples associated with a microarray experiment according to
#' the options in the config and combines them using the Robust Multichip
#' Average (RMA) algorithm, which background corrects, quantile normalizes, and
#' summarizes the data using median polish at the level of probes in an
#' expression set object (eset). If the option \code{normalization: no} is set
#' in the config, an expression set is returned which has neither been
#' background corrected nor normalized, a feature which is useful for quality
#' assessment.
#'
#' @param config Character string consisting of the path to the configuration
#' file generated using the \code{write.yaml.config} function or parsed
#' configuration list associated with a microarray experiment.
#'
#' @return Returns an expression set object that contains the microarray data at
#' the level of probes.
#'
#' @examples
#' ma.normalize(config)
#'
#' @references Irizarry, Rafael A., Bridget Hobbs, Francois Collin, Yasmin D.
#' Beazer-Barclay, Kristen J. Antonellis, Uwe Scherf, and Terence P. Speed.
#' "Exploration, normalization, and summaries of high density oligonucleotide
#' array probe level data." \emph{Biostatistics} 4, no. 2 (2003): 249-264.
#'
#' Carvalho, Benilton S., and Rafael A. Irizarry. "A framework for
#' oligonucleotide microarray preprocessing." \emph{Bioinformatics} 26, no. 19
#' (2010): 2363-2367.
#'
#' @seealso \code{\link{write.yaml.config}} to generate the configuration file,
#' and \code{\link{ma.summarize}} to convert an expression set object into a
#' list of genes and their associated log-fold changes and statistical values.
#'
#' @importFrom affyio read.celfile.header
#' @importFrom oligo read.celfiles rma
#' @importFrom parallel detectCores
#'
#' @export
ma.normalize <- function(config)
{
  config <- read.yaml.config(config)

  # Check that the cel files come from the same chip type
  check.cel.headers <- function(celFiles)
  {
    checkAtRandomPos <- sample(1:length(celFiles), min(10, length(celFiles)))
    randomSubsetChips <- vapply(celFiles[checkAtRandomPos], function(cel){return(affyio::read.celfile.header(cel)$cdfName)}, character(1))
    if(!all(randomSubsetChips[1] == randomSubsetChips))
    {
      subsetChips <- unique(randomSubsetChips)
      stop("Microarray samples come from different chip types.")
    }

    chipLibs <- .check.platform(randomSubsetChips[1])[c("platform.design","annotation.db")]

    return(chipLibs)
  }

  # Check that the platform design and annotation database packages are installed
  load.chip.libs <- function(chipLibs)
  {
    isInstalled <- TRUE
    isInstalled <- isInstalled & suppressMessages(requireNamespace(chipLibs[[1]], quietly = TRUE))
    isInstalled <- isInstalled & suppressMessages(requireNamespace(chipLibs[[2]], quietly = TRUE))
    if(!isInstalled)
    {
      stop(sprintf("Check that the packages '%s' and '%s' are installed from Bioconductor.", chipLibs[[1]], chipLibs[[2]]))
    }

    return(chipLibs)
  }

  # Read in all samples
  read.samples <- function(config)
  {
    cat("Reading in samples:")
    gfs <- .suppressOutput(oligo::read.celfiles, config$data$groups$sample.file)

    return(gfs)
  }

  # Normalize all samples using RMA
  norm.samples <- function(config, gfs)
  {
    normalize <- background <- ifelse(is.character(config$pipeline$normalization), TRUE, config$pipeline$normalization)
    filepath <- file.path(dirname(config$groups$group_file), ifelse(normalize, "pstEset.rds", "preEset.rds"))

    eset <- oligo::rma(gfs, background, normalize)

    # Save output to file
    if(config$global_options$save_intermediates == TRUE)
    {
      cat(sprintf("Saving ExpressionSet as '%s'... ", filepath))
      saveRDS(eset, filepath)
      cat("completed.\n")
    }

    return(eset)
  }

  # Enable multicore analysis
  Sys.setenv(R_THREADS = parallel::detectCores())

  # Suppress 'oligo' package diagnostics so that function output when reading samples isn't clobbered
  .silentLoad("oligo")

  load.chip.libs(check.cel.headers(config$data$groups$sample.file))

  rawGFS <- .gc.wrapper(read.samples, config = config)
  maEset <- .gc.wrapper(norm.samples, config = config, gfs = rawGFS)

  return(maEset)
}
