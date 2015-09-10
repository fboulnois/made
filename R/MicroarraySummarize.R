#' Microarray summarization
#'
#' Summarize a microarray experiment at the level of genes.
#'
#' Transforms an expression set (eset) which describes the microarray experiment
#' data at the probe level into a list of genes and their associated log-fold
#' changes and statistical values. The statistical values include the p-values,
#' the False Discovery Rate (FDR) adjusted p-values (q-values), and the 95\%
#' confidence intervals for the log-fold change. Expression sets which are not
#' log-transformed are log-transformed for the purpose of this function.
#' Batch effect is adjusted using Surrogate Variable Analysis (SVA) and
#' gene-level summarization is assessed using the empirical Bayes function from
#' the \code{limma} package.
#'
#' @param config Character string consisting of the path to the configuration
#' file generated using the \code{write.yaml.config} function or parsed
#' configuration list associated with a microarray experiment.
#'
#' @param eset Expression set object describing microarray experiment at the
#' level of probes.
#'
#' @return A list of top tables for each comparison of interest and their genes,
#' associated log-fold changes, and statistical values.
#'
#' @references Benjamini, Yoav, and Yosef Hochberg. "Controlling the false
#' discovery rate: a practical and powerful approach to multiple testing."
#' \emph{Journal of the Royal Statistical Society}. Series B (Methodological)
#' (1995): 289-300.
#'
#' Leek, Jeffrey T., and John D. Storey. "Capturing heterogeneity in gene
#' expression studies by surrogate variable analysis." \emph{PLoS Genet} 3, no.
#' 9 (2007): 1724-1735.
#'
#' Ritchie, Matthew E., Belinda Phipson, Di Wu, Yifang Hu, Charity W. Law, Wei
#' Shi, and Gordon K. Smyth. "limma powers differential expression analyses for
#' RNA-sequencing and microarray studies." \emph{Nucleic acids research} (2015):
#' gkv007.
#'
#' @importFrom AnnotationDbi keys select
#' @importFrom Biobase annotation exprs featureNames
#' @importFrom limma contrasts.fit eBayes lmFit topTable
#' @importFrom sva sva
#'
#' @export
ma.summarize <- function(config, eset)
{
  stopifnot(class(eset) == "ExpressionSet")
  config <- read.yaml.config(config)

  # Perform log transformation on eset if it has not already been done
  do.log2 <- function(eset)
  {
    curExprs <- Biobase::exprs(eset)
    qvals <- as.numeric(quantile(curExprs, c(0.00, 0.25, 0.50, 0.75, 0.99, 1.0), na.rm = TRUE))
    # Check whether or not data is log-transformed
    notLog <- (qvals[5] > 100) || (qvals[6] - qvals[1] > 50 && qvals[2] > 0) || (qvals[2] > 0 && qvals[2] < 1 && qvals[4] > 1 && qvals[4] < 2)
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
    hsel <- AnnotationDbi::select(hgdb, keys = AnnotationDbi::keys(hgdb), columns = c("PROBEID","ENTREZID","GENENAME","SYMBOL"))
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
      tt  <- limma::topTable(fit, coef = i, number = Inf, adjust.method = "BH", confint = TRUE)

      # Check which probe sets are differentially expressed or have significant log-fold change
      tt$DE <- tt$adj.P.Val < config$global_options$qvalue
      thresholdFC <- mean(tt$logFC) + 2*sd(tt$logFC)
      tt$sigFC <- abs(tt$logFC) > thresholdFC

      # Merge annotation database with differentially expressed probe sets, remove unknown probe sets
      tt$PROBEID <- row.names(tt)
      tt <- merge(tt, dbData, by = "PROBEID")
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
