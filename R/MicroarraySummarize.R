ma.summarize <- function(config, eset)
{
  stopifnot(class(eset) == "ExpressionSet")
  config <- read.yaml.config(config)

  # Perform log transformation on eset if it has not already been done
  do.log2 <- function(eset)
  {
    curExprs <- Biobase::exprs(eset)
    qvals <- as.numeric(quantile(curExprs, c(0.00, 0.25, 0.50, 0.75, 0.99, 1.0), na.rm=TRUE))
    # Check whether or not data is log-transformed
    notLog <- (qvals[5] > 100) || (qvals[6]-qvals[1] > 50 && qvals[2] > 0) || (qvals[2] > 0 && qvals[2] < 1 && qvals[4] > 1 && qvals[4] < 2)
    if(notLog)
    {
      warning("Expression set is not log-transformed.")
      curExprs[curExprs <= 0 ] <- NaN
      Biobase::exprs(eset) <- log2(curExprs)
    }
    return(eset)
  }

  # Load the corresponding annotation database
  get.annotation.data <- function(eset)
  {
    annotationPkg <- .check.platform(Biobase::annotation(eset))[["annotation.db"]]

    # Check that annotation package can be loaded
    isInstalled <- require(annotationPkg, character.only = TRUE, quietly = TRUE)
    if(!isInstalled)
    {
      stop(sprintf("Check that the package '%s' is installed from Bioconductor.", annotationPkg))
    }

    # Check that probes indicated in the eset are similar to those in the annotation database
    hgdb <- get(annotationPkg)
    hsel <- AnnotationDbi::select(hgdb, keys=AnnotationDbi::keys(hgdb), columns = c("PROBEID","ENTREZID","GENENAME","SYMBOL","UNIPROT"))
    if(length(setdiff(Biobase::featureNames(eset),hsel$PROBEID)) > 0)
    {
      stop("Expression set features must be probe sets.")
    }

    return(hsel)
  }

  # Adjust for surrogate variables (biological and non-biological variability)
  sva.estimate <- function(eset, groupData)
  {
    if(ncol(eset) != nrow(groupData$design.matrix))
    {
      stop("Number of samples in the expression set do not match those in the experimental design.")
    }

    modf <- groupData$design.matrix
    mod0 <- model.matrix(~1, groupData$groups)

    sva1 <- .gc.wrapper(sva::sva, Biobase::exprs(eset), modf, mod0)
    modSVs <- cbind(modf, sva1$sv)
    cat("\n")

    return(modSVs)
  }

  # Calculate differential expression for eset
  differential.expression <- function(eset, config, modSVs, dbData)
  {
    cmx <- config$data$contrast.matrix
    dmx <- config$data$design.matrix

    tf <- mf <- setNames(as.list(1:ncol(cmx)), colnames(cmx))
    for(i in 1:ncol(cmx))
    {
      cat(sprintf("Contrasting %s.\n", colnames(cmx)[i]))

      fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(eset, dmx), cmx))
      tt  <- limma::topTable(fit, coef = i, n = Inf, adjust.method = "BH", confint = TRUE)

      # Check which probe sets are differentially expressed or have significant log-fold change
      tt$DE <- tt$adj.P.Val < config$global_options$pvalue
      thresholdFC <- mean(tt$logFC) + 2*sd(tt$logFC)
      tt$sigFC <- abs(tt$logFC) > thresholdFC

      # Merge annotation database with differentially expressed probe sets, remove unknown probe sets
      tt$PROBEID <- row.names(tt)
      tt <- merge(tt, dbData, by="PROBEID")
      tt <- tt[!is.na(tt$ENTREZID), ]
      tt <- tt[ order(tt$P.Value) , ]

      # Store toptable and limma model in lists
      tf[[i]] <- tt
      mf[[i]] <- fit
    }

    return(list(top.table = tf, limma.model = mf, sva.model = modSVs))
  }

  eset <- do.log2(eset)

  dbData <- get.annotation.data(eset)
  modSVs <- sva.estimate(eset, config$data)
  models <- differential.expression(eset, config, modSVs, dbData)

  # Save output to file
  if(config$global_options$save_intermediates)
  {
    analysisDir <- dirname(config$groups$group_file)
    .write.table.list(models$top.table, analysisDir)
    saveRDS(models, file.path(analysisDir, "limmaModels.rds"))
  }

  return(models$top.table)
}
