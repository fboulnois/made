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
    validOpts <- c("pvalue", "save.intermediates",
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

    if(!is.null(opts$pvalue))
    {
      config$global_options$pvalue <- opts$pvalue
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

    return(oligoClasses::sampleNames(eset))
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
                              recursive = TRUE, ignore.case=TRUE)
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
