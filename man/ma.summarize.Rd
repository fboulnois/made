% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MicroarraySummarize.R
\name{ma.summarize}
\alias{ma.summarize}
\title{Microarray summarization}
\usage{
ma.summarize(config, eset)
}
\arguments{
\item{config}{Character string consisting of the path to the configuration
file generated using the \code{write.yaml.config} function or parsed
configuration list associated with a microarray experiment.}

\item{eset}{Expression set object describing microarray experiment at the
level of probes.}
}
\value{
A summary of all the data in the expression set for each group
comparison.
\describe{
   \item{top.tables}{A list of genes, associated log-fold changes, and
   other statistical values of interest.}
   \item{go.terms}{A list of biological terms associated with the expression
   set.}
   \item{reactome}{A list of biological pathways associated with the
   expression set.}
   \item{limma.model}{An \code{MArrayLM} fitted model object containing all
   of the statistical information relating to the expression set produced by
   the \code{eBayes} function in the \code{limma} package.}
   \item{design.matrix}{The experimental design matrix which represents the
   associations between samples and groups and is used to help fit the model.
   This matrix may have been modified by the \code{sva} function if batch
   effects were adjusted.}
 }
}
\description{
Summarize a microarray experiment at the level of genes.
}
\details{
Transforms an expression set (eset) which describes the microarray experiment
data at the probe level into a list of genes and their associated log-fold
changes and statistical values. The statistical values include the p-values,
the False Discovery Rate (FDR) adjusted p-values (q-values), and the 95\%
confidence intervals for the log-fold change. Expression sets which are not
log-transformed are log-transformed for the purpose of this function.
Batch effects are adjusted using Surrogate Variable Analysis (SVA) and
gene-level summarization is assessed using the empirical Bayes function from
the \code{limma} package. If the package \code{GOstats} is installed, Gene
Ontology (GO) term enrichment is performed to determine biologically relevant
terms in each group comparison. Similarly, if the \code{ReactomePA} package
is installed, pathway analysis is performed using the Reactome database to
determine which biological pathways are involved in each group comparison.
}
\examples{
if(require(madeData))
{
  dataPath <- system.file("extdata", package = "madeData")
  config <- file.path(dataPath, "config.yaml")
  eset <- readRDS(file.path(dataPath, "eset.rds"))
  ma.summarize(config, eset)
}

}
\references{
Benjamini, Yoav, and Yosef Hochberg. "Controlling the false
discovery rate: a practical and powerful approach to multiple testing."
\emph{Journal of the Royal Statistical Society}. Series B (Methodological)
(1995): 289-300.

Leek, Jeffrey T., and John D. Storey. "Capturing heterogeneity in gene
expression studies by surrogate variable analysis." \emph{PLoS Genet} 3, no.
9 (2007): 1724-1735.

Ritchie, Matthew E., Belinda Phipson, Di Wu, Yifang Hu, Charity W. Law, Wei
Shi, and Gordon K. Smyth. "limma powers differential expression analyses for
RNA-sequencing and microarray studies." \emph{Nucleic acids research} (2015):
gkv007.

Falcon, Seth, and Robert Gentleman. "Using GOstats to test gene lists for GO
term association." \emph{Bioinformatics} 23, no. 2 (2007): 257-258.

Croft, David, Gavin O'Kelly, Guanming Wu, Robin Haw, Marc Gillespie, Lisa
Matthews, Michael Caudy et al. "Reactome: a database of reactions, pathways
and biological processes." \emph{Nucleic acids research} (2010): gkq1018.
}

