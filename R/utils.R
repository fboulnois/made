# Try to load file if it exists otherwise perform function call
# Perform gc() after expensive load or function call
.gc.wrapper <- function(func, ...)
{
  on.exit(gc())
  return(do.call(func, list(...)))
}

# Verify that a microarray platform is supported
.check.platform <- function(platformName)
{
  platformData <- read.csv(system.file("extdata", "platforms.csv", package = "made"), stringsAsFactors = FALSE)

  pos <- which(platformData == platformName, arr.ind = TRUE)
  if(length(pos) == 0)
  {
    stop(sprintf("The microarray platform '%s' is not currently supported.", platformName))
  }

  return(platformData[pos[1], ])
}

# Versions of '*apply' functions with progress bar
.apply.internal <- function(applyFunc, X, FUN, ...)
{
  env <- environment()
  counter <- 0

  func.wrapper <- function(...)
  {
    curVal <- get("counter", envir = env)
    assign("counter", curVal + 1, envir = env)
    setTxtProgressBar(get("pbTxt", envir = env), curVal + 1)
    return(FUN(...))
  }

  argv <- list(X = X, FUN = func.wrapper, ...)
  if(as.character(applyFunc) == "apply")
  {
    if (!length(dim(X)))
    {
      stop("dim(X) must have a positive length")
    }
    pbLen <- ifelse(argv$MARGIN == 1, nrow(X), ncol(X))
  }
  else
  {
    pbLen <- length(X)
  }
  pbTxt <- txtProgressBar(min = 0, max = pbLen, style = 3)

  res <- do.call(deparse(applyFunc), argv)

  close(pbTxt)

  return(res)
}

.apply.pb <- function(X, MARGIN, FUN, ...)
{
  return(.apply.internal(quote(apply) , X = X, MARGIN = MARGIN, FUN = FUN, ...))
}

.lapply.pb <- function(X, FUN, ...)
{
  return(.apply.internal(quote(lapply), X = X, FUN = FUN, ...))
}

.sapply.pb <- function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE)
{
  return(.apply.internal(quote(sapply), X = X, FUN = FUN, ..., simplify = simplify, USE.NAMES = USE.NAMES))
}

.vapply.pb <- function(X, FUN, FUN.VALUE, ..., USE.NAMES = TRUE)
{
  return(.apply.internal(quote(vapply), X = X, FUN = FUN, FUN.VALUE = FUN.VALUE, ..., USE.NAMES = USE.NAMES))
}

# Split a vector into chunks
.chunk <- function(v, slices)
{
  return(split(v, ceiling(slices*seq_along(v)/length(v))))
}

# Full date timestamp
.timestamped <- function()
{
  return(paste0(gsub(":|-| ","", Sys.time())))
}

# Suppress output from a function
.suppressOutput <- function(func, ...)
{
  capture.output({
    res <- do.call(func, list(...))
  }, file = NULL)
  return(res)
}

# Write a list of tables into files
.write.table.list <- function(tableList, outputDir)
{
  stopifnot(is.list(tableList))
  stopifnot(is.character(outputDir))

  output.table <- function(tableList, tableNames, outputDir)
  {
    # Remove garbage characters
    tableFileName <- paste0(gsub("[^[:alnum:]]","", tableNames), ".csv")
    write.csv(tableList, file.path(outputDir, tableFileName), row.names = FALSE, na = "")
    return(NULL)
  }

  mapply(output.table, tableList = tableList, tableNames = names(tableList), outputDir = outputDir)
}
