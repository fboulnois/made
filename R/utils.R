# Perform gc() after expensive load or function call
.gc.wrapper <- function(func, ...)
{
  on.exit(gc())
  return(do.call(func, list(...)))
}

# Verify that a microarray platform is supported
.check.platform <- function(platformName)
{
  platformTrim <- tolower(gsub("[^[:alnum:]]","", platformName))
  platformData <- read.csv(system.file("extdata", "platforms.csv", package = "made"), stringsAsFactors = FALSE)

  pos <- which(platformData == platformName, arr.ind = TRUE)
  if(length(pos) == 0)
  {
    pos <- which(platformData == platformTrim, arr.ind = TRUE)
  }
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

# Rename columns of object and return modified object
.renameCols <- function(obj, oldCols, newCols)
{
  if(length(oldCols) == length(newCols))
  {
    colnames(obj)[na.omit(match(oldCols, colnames(obj)))] <- newCols
  }
  return(obj)
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
    tableFileName <- paste0(gsub("[^[:alnum:]-]","", tableNames), ".csv")
    write.csv(tableList, file.path(outputDir, tableFileName), row.names = FALSE, na = "")
    return(NULL)
  }

  mapply(output.table, tableList = tableList, tableNames = names(tableList), outputDir = outputDir)
}

# Load package and dependencies silently
.silentLoad <- function(package)
{
  loadPackage <- function(package)
  {
    # Code below based on 'library' function code
    libDirs <- .libPaths()
    libDirs <- libDirs[dir.exists(libDirs)]

    pkgPath <- find.package(package, libDirs, quiet = TRUE)
    if(length(pkgPath) == 0)
    {
      return(FALSE)
    }
    pkgDir  <- normalizePath(dirname(pkgPath), "/", TRUE)
    pkgInfo <- readRDS(system.file("Meta", "package.rds", package = package, lib.loc = pkgDir))

    .getRequiredPackages2(pkgInfo, quietly = TRUE)
    pkgDeps <- unique(names(pkgInfo$Depends))

    ns  <- loadNamespace(package, c(pkgDir, libDirs))
    env <- attachNamespace(ns, depends = pkgDeps)

    return(TRUE)
  }

  return(suppressMessages(tryCatch(loadPackage(package), error = function(e)
  {
    if(e$message == "namespace is already attached")
    {
      return(TRUE)
    }
    else
    {
      return(FALSE)
    }
  })))
}
