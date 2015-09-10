#' Write YAML configuration
#'
#' Generate the configuration file describing a microarray experiment.
#'
#' Generates the configuration and group file which describe a reproducible
#' microarray experiment. The files are created in the analysis directory using
#' the \code{analysisDir} option. The group file describes which samples belong
#' to which groups or conditions, for example normal controls and experimental
#' cases. The user must specify at which level they are grouping their samples
#' using the \code{groupBy} option:
#'
#'  \describe{
#'    \item{\code{"dirs" }}{Each directory containing samples can be assigned a
#'    group. This option is useful if there are many samples in the experiment.
#'    If this option is specified, before running this function the analysis
#'    directory must already contain subdirectories of samples.}
#'    \item{\code{"files"}}{Each individual sample can be assigned a group. If
#'    this option is specified, before running this function the analysis
#'    directory must already contain samples.}
#'    \item{\code{"eset" }}{Each sample specified in an expression set can be
#'    assigned a group. If this option is specified, the expression set must
#'    also be passed as an additional argument to the function.}
#'  }
#'
#' Other configuration options can also be specified as additional arguments,
#' otherwise defaults will be used:
#'
#'  \describe{
#'    \item{\code{qvalue}}{A number between 0 and 1 which corrects the pvalue
#'    for multiple comparisons and describes the minimum false discovery rate
#'    (FDR) at which a test may be called significant: the default is
#'    \code{0.05}.}
#'    \item{\code{save.intermediates}}{A logical indicating whether or not
#'    intermediate outputs should be saved: the default is \code{FALSE}.}
#'    \item{\code{quality.assessment}}{A logical indicating whether or not
#'    quality assessment should be performed: the default is \code{TRUE}.}
#'    \item{\code{normalization}}{Either \code{"rma"}, \code{TRUE}, or
#'    \code{FALSE}, indicating what sort of background correction and
#'    normalization should be performed, if any. Setting this value to
#'    \code{FALSE} is useful for quality assessment, since the expression set
#'    will not be background corrected or normalized: the default is
#'    \code{"rma"}.}
#'    \item{\code{summarization}}{A logical indicating whether or not gene level
#'    summarization should be performed: the default is \code{TRUE}.}
#'    \item{\code{write.report}}{A logical indicating whether or not a report of
#'    the experiment should be created: the default is \code{TRUE}.}
#'    \item{\code{group.file}}{A character string describing the name to give
#'    the group file: the default is to name the file \code{"groups.txt"}.}
#'    \item{\code{contrast.groups}}{A character string of one or more formulas
#'    separated by a comma describing how the different groups should be
#'    compared: by default groups \code{"A"}, \code{"B"}, and \code{"C"} are
#'    created and are each compared to each other.}
#'    \item{\code{eset}}{An expression set object describing the microarray
#'    experiment at the level of probes. This option should only be used if the
#'    \code{groupBy} option is set to \code{"eset"}, and thus has no default.}
#'    \item{\code{groups.df}}{A \code{data.frame} indicating which samples
#'    belong to which groups.
#'
#'    This option has no default: it can be used to programmatically set the
#'    groups for each sample ahead of time, otherwise the user will need to
#'    modify the group file after the configuration has been created.}
#'  }
#'
#' @param analysisDir Character string consisting of the path to the analysis
#' directory.
#'
#' @param groupBy Must be one of \code{"dirs"}, \code{"files"}, or
#' \code{"eset"}.
#'
#' @param ... Additional options to pass to the function.
#'
#' @seealso \code{\link{read.yaml.config}} to read and validate the
#' configuration.
write.yaml.config <- function(analysisDir, groupBy = NULL, ...)
{
  # Verify analysis directory and create it if it doesn't exist
  check.analysis.dir <- function(analysisDir)
  {
    if(!is.character(analysisDir))
    {
      stop("Parameter 'analysisDir' must be the path to a directory.")
    }
    if(!dir.exists(analysisDir))
    {
      dir.create(analysisDir, recursive = TRUE)
    }
  }

  # Verify that the 'groupBy' parameter is valid
  check.group.opt <- function(groupBy)
  {
    validParams <- c("dirs", "files", "eset")
    if(length(groupBy) != 1 || !is.element(groupBy, validParams))
    {
      txtParams <- paste0("'", validParams, "'", collapse = ", ")
      stop(sprintf("Parameter 'groupBy' must be one of %s.", txtParams))
    }
  }

  # Verify that the other specified options are valid and update the config if necessary
  check.other.opts <- function(groupBy, opts)
  {
    validOpts <- c("qvalue", "save.intermediates",
                   "quality.assessment", "normalization", "summarization", "write.report",
                   "group.file", "contrast.groups", "eset", "groups.df")

    pos <- is.element(names(opts), validOpts)
    if(!all(pos))
    {
      invalidOpt <- paste0("'", names(opts)[!pos][1], "'", collapse = ", ")
      stop(sprintf("Invalid option %s specified.", invalidOpt))
    }

    # Load default configuration file
    config <- normalizePath(system.file("extdata", "config.yaml", package = "made"), mustWork = TRUE)
    config <- yaml::yaml.load_file(config)

    if(!is.null(opts$qvalue))
    {
      config$global_options$qvalue <- opts$qvalue
    }

    if(!is.null(opts$save.intermediates))
    {
      config$global_options$save_intermediates <- opts$save.intermediates
    }

    if(!is.null(opts$quality.assessment))
    {
      config$pipeline$quality_assessment <- opts$quality.assessment
    }

    if(!is.null(opts$normalization))
    {
      config$pipeline$normalization <- opts$normalization
    }

    if(!is.null(opts$summarization))
    {
      config$pipeline$summarization <- opts$summarization
    }

    if(!is.null(opts$write.report))
    {
      config$pipeline$write_report <- opts$write.report
    }

    if(!is.null(opts$group.file))
    {
      config$groups$group_file <- opts$group.file
    }

    if(!is.null(opts$contrast.groups))
    {
      config$groups$contrast_groups <- opts$contrast.groups
    }

    config$groups$group_by <- groupBy

    # Check resulting configuration
    return(read.yaml.config(config, checkGroups = FALSE, getGroups = FALSE))
  }

  # Generate header for group file
  get.full.header <- function(groupType)
  {
    headerLines <- c(strwrap(
      sprintf(paste0("Directions: Assign groups, such as normal controls and experimental cases, ",
                     "to each sample %s by naming the elements in the 'Group' column. Names are ",
                     "alphanumeric strings and surrounded by double quotes. Each sample %s that ",
                     "belongs to the same group must share the same group name. By default, each ",
                     "sample %s is given the group name \"A\", \"B\", or \"C\", but any number of ",
                     "groups can be defined."),
              groupType, groupType, groupType),
      78, prefix = "# "), "# ",
      strwrap(paste0("Once groups are assigned, define how these groups will be compared in the ",
                     "file 'config.yaml' by modifying the 'contrast_groups' line. Each comparison ",
                     "should be separated by a comma and be composed of two groups separated by a ",
                     "dash."), 78, prefix = "# "))

    strSplitter <- "# --------------------------------------------------------------------------- #"
    return(paste(strSplitter, paste0(headerLines, collapse = "\n"), strSplitter, sep = "\n"))
  }

  # Check expression set and extract sample names from data
  get.eset.samples <- function(eset)
  {
    if(is.null(eset) || class(eset) != "ExpressionSet")
    {
      stop("Parameter groupBy = 'eset' was specified but a valid expression set was not passed to function.")
    }

    return(Biobase::sampleNames(eset))
  }

  # Based on whether groupBy = "dirs", "files", or "eset" take different actions
  get.group.data <- function(analysisDir, groupBy, opts)
  {
    if(groupBy == "dirs")
    {
      groupList <- list.dirs(analysisDir, full.names = FALSE)[-1]
      groupType <- "directory"
      middleCol <- "sample.dir"
    }
    else if(groupBy == "files")
    {
      groupList <- list.files(analysisDir, full.names = FALSE, pattern = "*.cel(.gz)?$",
                              recursive = TRUE, ignore.case = TRUE)
      groupType <- "file"
      middleCol <- "sample.file"
    }
    else if(groupBy == "eset")
    {
      groupList <- get.eset.samples(opts$eset)
      groupType <- "file"
      middleCol <- "sample.file"
    }

    if(length(groupList) == 0)
    {
      stop(sprintf("Parameter groupBy = %s was specified but no microarray %s were found in '%s'.",
                   groupBy, groupType, analysisDir))
    }

    if(!is.null(opts$groups.df))
    {
      df <- opts$groups.df
    }
    else
    {
      df <- setNames(data.frame(suppressWarnings(cbind(LETTERS[1:min(length(groupList),3)], groupList, basename(groupList)))),
                     c("Group", middleCol, "description"))
    }

    return(list(header = get.full.header(groupType), df = df))
  }

  # Write config file output
  write.config.file <- function(analysisDir, config)
  {
    configFile <- file.path(analysisDir, "config.yaml")
    writeLines(yaml::as.yaml(config), con = configFile)
  }

  # Write group file output
  write.group.file <- function(analysisDir, config, group)
  {
    groupFile <- file.path(analysisDir, config$groups$group_file)
    groupHndl <- file(groupFile, open = "wt")
    on.exit(close(groupHndl))

    writeLines(group$header, con = groupHndl)
    write.table(group$df, groupHndl, row.names = FALSE, sep = "\t")

    # Give user warning of what to do next
    warning(sprintf("The file '%s' must be modified to assign controls and cases before launching analysis!",
                    normalizePath(groupFile)), immediate. = TRUE)
  }

  check.analysis.dir(analysisDir)
  check.group.opt(groupBy)

  opts <- list(...)

  config <- check.other.opts(groupBy, opts)
  groups <- get.group.data(analysisDir, groupBy, opts)

  write.config.file(analysisDir, config)
  write.group.file(analysisDir, config, groups)
}
