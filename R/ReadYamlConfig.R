read.yaml.config <- function(config, checkGroups = TRUE, getGroups = TRUE)
{
  # If YAML has already been validated avoid further processing
  if(class(config) == "MadeConfig")
  {
    return(config)
  }

  # Try to read file or string as YAML
  cfgdir <- NULL
  if(is.character(config))
  {
    if(file.exists(config))
    {
      cfgdir <- dirname(config)
      config <- yaml::yaml.load_file(config)
    }
    else
    {
      config <- yaml::yaml.load(config)
    }
  }

  # Error out if config couldn't be read as YAML
  if(!is.list(config))
  {
    stop("Parameter 'config' must be a YAML file, YAML string, or YAML object containing options for microarray analysis.")
  }

  # Get the original name of a parameter
  get.varname <- function(param)
  {
    for(n in 1:sys.nframe())
    {
      varName <- deparse(substitute(param, env = parent.frame(n)))
      if(varName != "param")
      {
        break
      }
    }
    return(varName)
  }

  # If condition is not met, error out with message
  cerr <- function(cond, param, txt = NULL, error = TRUE)
  {
    if(!cond && error)
    {
      stop(ifelse(is.null(txt), sprintf("Missing configuration section '%s'.", get.varname(param)),
                  sprintf("Expected configuration section '%s' to be %s, but got '%s' instead.",
                          get.varname(param), txt, param)))
    }
    return(cond)
  }

  # Functions to check that YAML nodes match certain conditions
  not.null <- function(param, error = TRUE)
  {
    return(cerr(!is.null(param), param = param, error = error))
  }

  is.single.logical <- function(param, error = TRUE)
  {
    return(cerr(not.null(param) && length(param) == 1 && is.logical(param), param, "either yes or no", error))
  }

  is.probability <- function(param, error = TRUE)
  {
    return(cerr(not.null(param) && length(param) == 1 &&
                is.double(param) && param >= 0 && param <= 1, param, "a number between 0 and 1", error))
  }

  is.single.string <- function(param, error = TRUE)
  {
    return(cerr(not.null(param) && length(param) == 1 && is.character(param), param, "a string", error))
  }

  is.file <- function(param, filename, error = TRUE)
  {
    return(cerr(is.single.string(filename) && file.exists(filename), param, "a valid file", error))
  }

  is.option <- function(param, opts, error = TRUE)
  {
    txt <- sprintf("one of %s", paste0("'", opts, "'", collapse = ", "))
    return(cerr(is.single.string(param) && is.element(param, opts), param, txt, error))
  }

  is.valid.norm <- function(param)
  {
    cerr(is.single.logical(param, error = FALSE) || is.option(param, "rma", error = FALSE),
         param, "one of 'rma', yes, or no")
  }

  # Verify that global options section is valid
  check.goptions <- function(config)
  {
    not.null(config$global_options)
    is.probability(config$global_options$pvalue)
    is.single.logical(config$global_options$save_intermediates)
  }

  # Verify that pipeline section is valid
  check.pipeline <- function(config)
  {
    not.null(config$pipeline)
    is.single.logical(config$pipeline$quality_assessment)
    is.valid.norm(config$pipeline$normalization)
    is.single.logical(config$pipeline$summarization)
    is.single.logical(config$pipeline$write_report)
  }

  # Verify that groups section is valid
  check.groups <- function(config)
  {
    not.null(config$groups)
    is.option(config$groups$group_by, c("dirs", "files", "eset"))
    is.single.string(config$groups$contrast_groups)
  }

  # Special validation for group file
  check.group.file <- function(param, cfgdir, checkGroups)
  {
    # Get full path if local path is specified
    fileName <- ifelse(dirname(param) == "." && !is.null(cfgdir), file.path(cfgdir, param), param)
    is.file(param, fileName, checkGroups)
    return(ifelse(checkGroups, normalizePath(fileName, mustWork = FALSE), fileName))
  }

  # Check configuration sections
  check.goptions(config)
  check.pipeline(config)
  check.groups(config)

  config$groups$group_file <- check.group.file(config$groups$group_file, cfgdir, checkGroups)

  if(getGroups)
  {
    config <- read.group.file(config)
    class(config) <- "MadeConfig"
  }

  return(config)
}