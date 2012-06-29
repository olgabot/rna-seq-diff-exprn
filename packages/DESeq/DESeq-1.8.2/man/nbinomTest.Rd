\name{nbinomTest}
\Rdversion{1.1}
\alias{nbinomTest}

\title{
   Test for differences between the base means for two conditions
}
\description{
   This function tests for differences between the base means of two
   conditions (i.e., for differential expression in the case of 
   RNA-Seq). 
}
\usage{
nbinomTest(cds, condA, condB, pvals_only = FALSE, eps=NULL )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cds}{
      a CountDataSet with size factors and raw variance functions
}
  \item{condA}{
      one of the conditions in 'cds'
}
  \item{condB}{
      another one of the conditions in 'cds'
}
  \item{pvals_only}{
      return only a vector of (unadjusted) p values instead of the data frame described below.
}
  \item{eps}{
     This argument is no longer used. Do not use it.
}
}
\details{
   See \code{\link{nbinomTestForMatrices}} for more technical informations
}
\value{
A data frame with the following columns:
\item{id}{The ID of the observable, taken from the row names of the counts slots.}
\item{baseMean}{The base mean (i.e., mean of the counts divided by the size factors) for the counts 
for both conditions}
\item{baseMeanA}{The base mean (i.e., mean of the counts divided by the size factors) for the counts 
for condition A}
\item{baseMeanB}{The base mean for condition B}
\item{foldChange}{The ratio meanB/meanA}
\item{log2FoldChange}{The log2 of the fold change}
\item{pval}{The p value for rejecting the null hypothesis 'meanA==meanB'}
\item{padj}{The adjusted p values (adjusted with 'p.adjust( pval, method="BH")')}

}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds )
head( nbinomTest( cds, "A", "B" ) )
}
