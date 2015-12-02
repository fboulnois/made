#' Read group file
#'
#' Read and validate the group file associated with a microarray experiment.
#'
#' \strong{Note}: This function is mostly for internal use. Users should use
#' \code{read.yaml.config} instead because it calls this function by default.
#'
#' The group file describes which samples belong to which groups or conditions,
#' for example normal controls and experimental cases. Groups must be assigned
#' prior to launching the analysis. Groups are compared according to the option
#' \code{contrast_groups} in the configuration and can be defined at three
#' different levels:
#'
#'  \describe{
#'    \item{\code{"dirs" }}{Each directory containing samples can be assigned a
#'    group.}
#'    \item{\code{"files"}}{Each individual sample can be assigned a group.}
#'    \item{\code{"eset" }}{Each sample specified in an expression set can be
#'    assigned a group.}
#'  }
#'
#' @param config Character string consisting of the path to the configuration
#' file generated using the \code{write.yaml.config} function or parsed
#' configuration list associated with a microarray experiment.
#'
#' @return Returns a modified configuration with group data.
#'
#' @examples
#' if(require(madeData))
#' {
#'   config <- system.file("extdata", "config.yaml", package = "madeData")
#'   read.group.file(config)
#' }
#'
#' @seealso \code{\link{write.yaml.config}} to generate the configuration file
#' and \code{\link{read.yaml.config}} to read and validate it.
#'
#' @importFrom limma makeContrasts
#' @importFrom stringr str_extract_all str_replace_all str_split str_trim
#'
#' @export
read.group.file <- function(config)
{
  config <- read.yaml.config(config, getGroups = FALSE)

  # Give error when invalid groups are specified in group file
  invalid.groups.error <- function(df, groupType, invalidGroups)
  {
    invalidGroupsTxt <- paste0("'", invalidGroups, "'", collapse = ", ")
    stop(sprintf("The %s %s specified in the file '%s' do not exist in the directory '%s'.",
                 groupType, invalidGroupsTxt, basename(attr(df, ".filename")), dirname(attr(df, ".filename"))))
  }

  # Verify that the group file is properly formatted
  check.group.file <- function(df, groupType, middleCol, callFunc, fileCheck = TRUE)
  {
    if(!identical(colnames(df), c("Group", middleCol, "description")))
    {
      stop(sprintf("The file '%s' must contain a tab-seperated matrix with the column names \"Group\", \"%s\", and \"description\". Lines starting with '#' are ignored.", basename(attr(df, ".filename")), middleCol))
    }
    if(nrow(df) < 2)
    {
      stop(sprintf("The file '%s' must have at least 2 or more rows.", basename(attr(df, ".filename"))))
    }
    if(!is.character(df$Group))
    {
      stop(sprintf("The data of the 'Group' column in the file '%s' must all be alphanumeric strings.",
                   basename(attr(df, ".filename"))))
    }
    if(!is.character(df[[middleCol]]))
    {
      stop(sprintf("The data of the '%s' column in the file '%s' must all be %s which exist in the directory '%s'.",
                   middleCol, basename(attr(df, ".filename")), groupType, dirname(attr(df, ".filename"))))
    }
    groups <- file.path(dirname(attr(df, ".filename")), df[[middleCol]])

    groupsExist <- vapply(groups, callFunc, logical(1))
    if(fileCheck && !all(groupsExist))
    {
      invalid.groups.error(df, groupType, basename(groups[!groupsExist]))
    }
  }

  # Get the full paths and groups of each of the files from their directories
  get.dir.groups <- function(df)
  {
    celFiles <- list.files(path = dirname(attr(df, ".filename")), pattern = "*.cel(.gz)?$",
                           full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    parentDirs <- basename(dirname(celFiles))

    pos <- match(parentDirs, df$sample.dir)

    new.df <- data.frame(Group = df$Group[pos], sample.file = normalizePath(celFiles),
                         description = df$description[pos], stringsAsFactors = FALSE)

    attr(new.df, ".filename") <- attr(df, ".filename")

    return(new.df)
  }

  # Get the full paths of each of the files
  get.file.groups <- function(df)
  {
    df$sample.file <- normalizePath(file.path(dirname(attr(df,".filename")), df$sample.file))
    return(df)
  }

  # Verify that the groups exist and are well formatted
  check.groups <- function(df, contrastLine)
  {
    cnlTokens <- stringr::str_extract_all(contrastLine, "\\w+")[[1]]
    badTokens <- setdiff(cnlTokens, df$Group)
    if(length(badTokens) != 0)
    {
      badTokensTxt <- paste0("'", badTokens, "'", collapse = ",")
      stop(sprintf("Contrast groups %s defined in 'config.yaml' do not exist in '%s'", badTokensTxt, basename(attr(df, ".filename"))))
    }

    cnlFormulas <- stringr::str_trim(stringr::str_split(contrastLine, ",")[[1]])

    dmx <- model.matrix(~0 + Group, df)
    colnames(dmx) <- stringr::str_replace_all(colnames(dmx), "Group", "")

    cmx <- tryCatch(
      {
        do.call(limma::makeContrasts, c(as.list(cnlFormulas), list(levels = as.factor(df$Group))))
      }, error = function(e)
      {
        badFormulasTxt <- paste0("'", cnlFormulas, "'", collapse = ", ")
        stop(sprintf("One or more of the following contrast groups in the file 'config.yaml' are not valid: %s.", badFormulasTxt))
      })

    return(list(design.matrix = dmx, contrast.matrix = cmx))
  }

  # Based on whether groupBy == "dirs", "files", "eset" take different actions
  read.group.table <- function(config)
  {
    df <- read.table(config$groups$group_file, sep = "\t", header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
    attr(df, ".filename") <- config$groups$group_file

    if(config$groups$group_by == "dirs")
    {
      groupType <- "directories"
      middleCol <- "sample.dir"
      callFunc  <- quote(dir.exists)

      check.group.file(df, groupType, middleCol, callFunc)

      df <- get.dir.groups(df)
    }
    else if(config$groups$group_by == "files")
    {
      groupType <- "files"
      middleCol <- "sample.file"
      callFunc  <- quote(file.exists)

      check.group.file(df, groupType, middleCol, callFunc)

      df <- get.file.groups(df)
    }
    else if(config$groups$group_by == "eset")
    {
      groupType <- "files"
      middleCol <- "sample.file"
      callFunc  <- quote(file.exists)

      check.group.file(df, groupType, middleCol, callFunc, fileCheck = FALSE)
    }

    expMatrices <- check.groups(df, config$groups$contrast_groups)

    return(list(groups = df, design.matrix = expMatrices$design.matrix, contrast.matrix = expMatrices$contrast.matrix))
  }

  config$data <- read.group.table(config)

  return(config)
}
