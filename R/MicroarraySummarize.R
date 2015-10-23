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
  sva.estimate <- function(eset, config)
  {
    modf <- config$data$design.matrix
    mod0 <- model.matrix(~1, config$data$groups)

    cmx <- config$data$contrast.matrix

    if(ncol(eset) != nrow(modf))
    {
      stop("Number of samples in the expression set do not match those in the experimental design.")
    }

    if(config$global_options$adjust_batch_effect)
    {
      sva1 <- .gc.wrapper(sva::sva, Biobase::exprs(eset), modf, mod0)
      cat("\n")

      modf <- cbind(modf, sva1$sv)

      dummy <- matrix(0, nrow = sva1$n.sv, ncol = ncol(cmx))
      cmx <- rbind(cmx, dummy)
    }

    return(list(design.matrix = modf, contrast.matrix = cmx))
  }

  # Calculate differential expression for eset
  differential.expression <- function(eset, config, models, dbData)
  {
    dmx <- models$design.matrix
    cmx <- models$contrast.matrix

    # Compute microarray statistics using empirical Bayes
    fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(eset, dmx), cmx))
    ess <- fit$coefficients / sqrt(fit$s2.post)

    hasGODB <- .silentLoad("GO.db")
    hasKEGG <- .silentLoad("KEGGREST")

    tf <- goid <- kegg <- setNames(vector(mode = "list", length = ncol(cmx)), colnames(cmx))
    for(i in 1:ncol(cmx))
    {
      cat(sprintf("Contrasting %s.\n", colnames(cmx)[i]))

      tt <- limma::topTable(fit, coef = i, number = Inf, adjust.method = "BH", confint = TRUE)
      tt$PROBEID <- row.names(tt)

      # Merge effect size and variance for each probe
      tt <- merge(tt, data.frame(PROBEID = names(fit$s2.post), ES = ess[ ,i], s2 = fit$s2.post), by = "PROBEID")

      # Check which probe sets are differentially expressed or have large log-fold change
      tt$DE <- tt$adj.P.Val < config$global_options$qvalue
      tt$sigFC <- abs(tt$logFC) > log2(1.5)

      # Merge annotation database with differentially expressed probe sets, remove unknown probe sets
      tt <- merge(tt, dbData, by = "PROBEID")
      tt <- tt[!is.na(tt$ENTREZID), ]
      tt <- tt[ order(tt$P.Value) , ]
      tt <- tt[!duplicated(tt$ENTREZID), ]

      # Store toptable, go terms, and kegg pathways in lists
      tf[[i]] <- tt

      if(hasGODB)
      {
        goidRes <- limma::goana(tt$ENTREZID[tt$DE & tt$sigFC], universe = tt$ENTREZID, covariate = tt$AveExpr)
        goid[[i]] <- goidRes[order(goidRes$P.DE), ]
      }
      if(hasKEGG)
      {
        keggRes <- limma::kegga(tt$ENTREZID[tt$DE & tt$sigFC], universe = tt$ENTREZID, covariate = tt$AveExpr)
        kegg[[i]] <- keggRes[order(keggRes$P.DE), ]
      }
    }

    return(list(top.tables = tf, go.terms = goid, kegg.pathways = kegg, limma.model = fit, design.matrix = dmx))
  }

  eset <- .do.log2(eset)

  dbData <- get.annotation.data(eset)
  models <- sva.estimate(eset, config)
  models <- differential.expression(eset, config, models, dbData)

  # Save output to file
  if(config$global_options$save_intermediates)
  {
    analysisDir <- dirname(config$groups$group_file)
    .write.table.list(models$top.tables, analysisDir)
    saveRDS(models, file.path(analysisDir, "limmaModels.rds"))
  }

  return(models$top.tables)
}
