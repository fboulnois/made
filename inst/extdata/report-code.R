# ---- init ----

figNum <- 0
draw.fig <- function(func, ...)
{
  do.call(func, list(...))
  figNum <<- figNum + 1
  return(figNum)
}

hasGoTerms <- !is.null(results$go.terms)
hasReactPA <- !is.null(results$reactome)

# ---- optiontext ----

logi.to.str <- function(opt, plural = FALSE)
{
  txt <- ""
  if(opt)
  {
    if(plural)
    {
      txt <- "**were**"
    }
    else
    {
      txt <- "**was**"
    }
  }
  else
  {
    if(plural)
    {
      txt <- "were **not**"
    }
    else
    {
      txt <- "was **not**"
    }
  }
  return(txt)
}

# ---- sampleinfo ----

sample.table <- function(config)
{
  group.samples <- function(group, groupData)
  {
    pos <- group == groupData$Group
    txt <- paste0(groupData[pos,"sample.file"], collapse = ", ")
    return(txt)
  }
  groupData <- config$data$groups
  uniGroups <- unique(groupData$Group)
  res <- vapply(uniGroups, group.samples, character(1), groupData = groupData)
  df <- data.frame(Group = uniGroups, Samples = res)
  print(knitr::kable(df, row.names = FALSE))
  cat("\n\n")
  return(list(groups = uniGroups, samples = groupData$sample.file))
}

sampleInfo <- sample.table(config)

# ---- plothistogram ----

date.parse <- function(dates)
{
  # Check what date format a column could be
  col.format <- function(col)
  {
    vals <- as.numeric(col)
    fmts <- c("%y", "%m", "%d")
    if(any(vals > 31))
    {
      fmts[[2]] <- NA
      fmts[[3]] <- NA
    }
    else if(any(vals > 12))
    {
      fmts[[2]] <- NA
    }
    bads <- is.na(fmts)
    return(fmts[!bads])
  }

  if(is.null(dates))
  {
    stop("Could not read in date vector")
  }

  # Try to extract general date format
  tknDat <- stringr::str_match(dates, "(\\d+)[-/](\\d+)[-/](\\d+)")
  tknSep <- stringr::str_extract(dates[1],"[-/]")
  if(any(is.na(tknDat)))
  {
    stop("Could not determine date string format")
  }

  # Determine format and frequency of each date column
  dcol <- apply(tknDat[,2:4], 2, col.format)
  freq <- apply(tknDat[,2:4], 2, function(col){ return(length(unique(col))) } )

  # Second position should never be year
  dcol[[2]] <- setdiff(dcol[[2]], "%y")

  # Most frequent position should be day, least frequent should be year
  dpos <- freq == max(freq)
  if(sum(dpos) == 1)
  {
    dcol[[which(dpos)]] <- "%d"
  }
  ypos <- freq == min(freq)
  if(sum(ypos) == 1)
  {
    dcol[[which(ypos)]] <- "%y"
  }

  # Each date column format that is certain is removed from the other columns
  for(i in 1:3)
  {
    zcol <- dcol[[i]]
    if(length(zcol) == 1)
    {
      pos1 <- (i + 0) %% 3 + 1
      pos2 <- (i + 1) %% 3 + 1

      dcol[[pos1]] <- setdiff(dcol[[pos1]], zcol)
      dcol[[pos2]] <- setdiff(dcol[[pos2]], zcol)
    }
  }

  # Error if after format removal there are still columns with more than one format
  notLen1 <- vapply(dcol, length, numeric(1)) != 1
  if(any(notLen1))
  {
    stop("Could not resolve date/column relationship")
  }
  else
  {
    fmtVec <- unlist(dcol)
  }

  # Check whether the date provides a 4-digit year
  ypos <- which(fmtVec == "%y")
  if(nchar(tknDat[1,ypos + 1]) == 4)
  {
    fmtVec[ypos] <- "%Y"
  }

  # Concatenate various date format pieces together
  fmtStr <- gsub("[-/]$", "", paste0(fmtVec, tknSep, collapse = ""))

  return(as.Date(tknDat[,1], fmtStr))
}

# Take the rolling mean of consecutive numbers
rollmean <- function(centers, i)
{
  centers <- as.vector(centers)
  if(i == 1)
  {
    return(centers)
  }
  else
  {
    n <- length(centers) - 1
    v <- rep(0, n)
    for(j in 1:n)
    {
      v[j] <- (centers[j] + centers[j + 1])/2
    }
    return(v)
  }
}

# Determine how many groups of expression set scan dates exist
get.batches <- function(eset)
{
  scanDates <- date.parse(Biobase::pData(Biobase::protocolData(eset))$dates)
  minDate <- min(scanDates)

  days <- as.numeric(scanDates - minDate)

  for(i in 2:10)
  {
    kmn <- kmeans(days, centers = i, nstart = 3)
    pct <- kmn$betweenss / kmn$totss
    sep <- c(0, rollmean(kmn$centers, i), max(days) + 1)

    if(pct > 0.90)
    {
      break
    }
  }

  lvls <- levels(cut(days, sep, include.lowest = TRUE, dig.lab = 4))
  text <- paste(gsub("\\.\\d+", "", lvls), "days")

  return(list(n = i, text = text, kmeans = kmn))
}

draw.histogram <- function(eset)
{
  # Try to extract batches by scan date
  hasBatches <- FALSE
  batchColors <- tryCatch({
    batches <- get.batches(eset)
    hasBatches <- TRUE
    return(cm.colors(batches$n)[batches$kmeans$cluster])
  }, error = function(e)
  {
    sampleNum <- length(Biobase::sampleNames(eset))
    return(cm.colors(sampleNum))
  })

  # Draw the histogram
  draw.fig(oligo::hist, x = eset, col = batchColors, main = "Log-intensity values distribution")
  if(hasBatches)
  {
    legend("topright", legend = batches$text, fill = cm.colors(batches$n), inset = 0.01, bg = "white")
  }

  return(hasBatches)
}

hasBatches <- draw.histogram(eset)

# ---- plotdendro  ----

color.leaf <- function(node, samples, sampleColors)
{
  if(is.leaf(node))
  {
    gsm <- stringr::str_extract(attr(node, "label"), "[\\w-]+")
    pos <- which(gsm == samples)
    attr(node, "nodePar") <- list(lab.col = sampleColors[pos])
  }
  return(node)
}

draw.dendrogram <- function(config, eset)
{
  # Color samples by their groups
  groups <- unique(config$data$groups$Group)
  posGroups <- match(config$data$groups$Group, groups)
  dendroColors <- rainbow(length(groups))

  # Load expression matrix and replace na values
  ex <- Biobase::exprs(eset)
  ex[is.na(ex)] <- 0

  # Spearman correlation dendrogram of sample data
  dS <- as.dist(1 - cor(ex, method = "spearman"))
  sampleClustering <- hclust(dS)
  sampleDendrogram <- as.dendrogram(sampleClustering, 0.05)
  sampleDendrogram <- dendrapply(sampleDendrogram, color.leaf, samples = config$data$groups$sample.file, sampleColors = dendroColors[posGroups])

  par(cex = 0.8)
  draw.fig(plot, x = sampleDendrogram, main = "Hierarchical clustering of samples")
  legend("topright", legend = groups, fill = dendroColors, inset = 0.01, bg = "white")
}

draw.dendrogram(config, eset)

# ---- degenestbl ----

top50.genes <- function(results)
{
  allNames <- names(results$top.tables)
  numNames <- length(allNames)
  for(i in 1:numNames)
  {
    maxRows <- 1:min(nrow(results$top.tables[[i]]), 50)
    tt <- results$top.tables[[i]][maxRows, c("PROBEID", "ENTREZID", "SYMBOL", "logFC", "CI.L", "CI.R", "adj.P.Val")]
    tt$qvalue <- ifelse(tt$adj.P.Val < 1e-16, "< 1e-16", format(signif(tt$adj.P.Val, 3), scientific = TRUE))
    tt$CI <- sprintf("%.3f to %.3f", tt$CI.L, tt$CI.R)
    tt$SYMBOL <- sprintf("[%s](http://www.ncbi.nlm.nih.gov/gene/%s)", tt$SYMBOL, tt$ENTREZID)
    rownames(tt) <- NULL

    tt <- tt[, c("SYMBOL", "logFC", "CI", "qvalue")]
    colnames(tt) <- c("Gene symbol", "Log fold-change", "95% CI", "qvalue")

    cat(sprintf("## %s\n", allNames[[i]]))
    print(knitr::kable(tt, digits = 3, align = c("l", "r", "r", "r")))
    cat("\n\n")
  }
}

top50.genes(results)

# ---- corheatmap ----

cor.heatmap <- function(eset, results)
{
  allNames <- names(results$top.tables)
  numNames <- length(allNames)
  for(i in 1:numNames)
  {
    tt <- results$top.tables[[i]]

    # Calculate correlations of most differentially expressed probes
    numProbes <- min(2000, length(tt$PROBEID))
    topProbes <- tt$PROBEID[1:numProbes]
    dS <- cor(t(Biobase::exprs(eset[topProbes, ])))

    # Find most correlated genes from probes
    numGenes <- min(ncol(dS), 50)
    mainTitle <- sprintf("Correlation heatmap of the %d\nmost differentially expressed genes in %s", numGenes, allNames[[i]])
    corPos <- !is.na(dS) & dS > 0.9 & dS < 0.9999
    corTop <- sort(colSums(corPos), decreasing = TRUE)
    corProbes <- names(corTop[1:numGenes])
    corGenes <- tt$SYMBOL[tt$PROBEID %in% corProbes]

    # Draw correlation heatmap of genes
    hmColors <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
    cat(sprintf("## %s\n", allNames[[i]]))
    par(oma = c(0, 0, 2.2, 0))
    draw.fig(heatmap, x = dS[corProbes, corProbes], Rowv = NA, Colv = NA, symm = TRUE, labRow = corGenes, labCol = corGenes, main = mainTitle, xlab = "genes", ylab = "genes", col = hmColors, margins = c(7,7))
    cat("\n\n")
  }
}

cor.heatmap(eset, results)

# ---- goterms ----

if(hasGoTerms)
{
  ontogo <- list(BP = "Biological Process", CC = "Cellular Component", MF = "Molecular Function")

  print.goterms <- function(goresult)
  {
    goterms <- GOstats::summary(goresult)
    if(nrow(goterms) == 0)
    {
      return(NULL)
    }
    maxRows <- 1:min(nrow(goterms), 20)
    goidcol <- colnames(goterms)[[1]]
    goterms <- goterms[maxRows, c(goidcol, "Term", "OddsRatio", "Pvalue")]

    goterms$Pvalue <- ifelse(goterms$Pvalue < 1e-16, "< 1e-16", format(signif(goterms$Pvalue, 3), scientific = TRUE))
    goterms$OddsRatio <- ifelse(goterms$OddsRatio > 500, "> 500.00", format(goterms$OddsRatio, digits = 3))
    goterms[[goidcol]] <- sprintf("[%s](http://amigo.geneontology.org/amigo/term/%s)", goterms[[goidcol]], goterms[[goidcol]])

    colnames(goterms) <- c("GO ID", "Term", "Odds Ratio", "Conditional pvalue")

    ontoWhich <- stringr::str_match(goidcol, "GO(\\w{2})ID")[[2]]
    ontoPos <- match(ontoWhich, names(ontogo))

    cat(sprintf("#### %s\n", ontogo[ontoPos]))
    print(knitr::kable(goterms, digits = 3, align = c("l", "l", "r", "r")))
    cat("\n\n")
  }

  top20.goterms <- function(results)
  {
    allNames <- names(results$go.terms)
    numNames <- length(allNames)
    for(i in 1:numNames)
    {
      if(is.null(results$go.terms[[i]]) || (is.atomic(results$go.terms[[i]]) && is.na(results$go.terms[[i]])))
      {
        next
      }
      cat(sprintf("### %s\n", allNames[[i]]))
      lapply(results$go.terms[[i]], print.goterms)
      cat("\n\n")
    }
  }

  top20.goterms(results)
}

# ---- reactpa ----

if(hasReactPA)
{
  top20.reactome <- function(results)
  {
    allNames <- names(results$reactome)
    numNames <- length(allNames)
    for(i in 1:numNames)
    {
      reactome <- DOSE::summary(results$reactome[[i]])
      if(is.null(results$reactome[[i]]) || (is.atomic(results$reactome[[i]]) && is.na(results$reactome[[i]])) || nrow(reactome) == 0)
      {
        next
      }
      maxRows <- 1:min(nrow(reactome), 20)
      reactome <- reactome[maxRows, c("ID", "Description", "qvalue")]
      reactome$qvalue <- ifelse(reactome$qvalue < 1e-16, "< 1e-16", format(signif(reactome$qvalue, 3), scientific = TRUE))
      reactome$ID <- sprintf("[%s](http://www.reactome.org/content/detail/%s)", reactome$ID, reactome$ID)
      rownames(reactome) <- NULL

      cat(sprintf("### %s\n", allNames[[i]]))
      print(knitr::kable(reactome, digits = 3, align = c("l", "l", "r")))
      cat("\n\n")
    }
  }

  topReactome <- top20.reactome(results)
}

# ---- sessioninfo ----

sessionInfo()
