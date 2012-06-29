\name{newCountDataSet}
\Rdversion{1.1}
\alias{newCountDataSet}

\title{
Create a CountDataSet object
}
\description{
This function creates a CountDataSet object from a matrix or data frame of count data.
}
\usage{
newCountDataSet(countData, conditions, sizeFactors = NULL, phenoData = NULL, featureData = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{countData}{
      A matrix or data frame of count data, i.e., of non-negative integer values. The rows correspond
      to observations (e.g., number of reads that were assigned to a gene), the columns correspond
      to samples (or experiments). Note that biological replicates should each get their own column,
      while the counts of technical replicates (i.e., several sequencing ruins/lanes from the same
      sample) have to be summed up into a single column.
}
  \item{conditions}{
      A factor of experimental conditions (or treatments, or tissue types, or phenotypes, or the like).
      The length of the factor has to be equal to the number of columns of the countData matrix,
      assigning a condition to each sample. If 'conditions' is not a factor, it will be converted to one.
      
      Alternatively, you may pass a data frame, that will be placed in pData(cds) as is and can then
      be used with the modes "pooled" and "blind" in \code{\link{estimateVarianceFunctions}} and its
      columns can be refered top in a model formula provided to \code{\link{fitNbinomGLMs}}.
}
  \item{sizeFactors}{
      This argument is deprecated. Do not use it. (Size factors should always be estimated from
      the data with \code{\link{estimateSizeFactors}}. If you need to set size factors manually for
      some reasons, change the \code{pData(cds)$sizeFactor}.
}
  \item{phenoData}{
      You may pass an AnnotatedDataFrame here to describe the columns of the count matrix. Note that the
      package always adds two rows (or creates a new AnnotatedDataFrame with only these two rows in case
      you do not supply one) with names "condition" and "sizeFactor" to store this information.
}
  \item{featureData}{
      You may pass an AnnotatedDataFrame here to describe the rows of the count matrix. The package will just
      pass through this information without using it. Note that further columns will be added to feature data
      later, when estimating dispersions.
}
}
\details{
   See also \code{\link{CountDataSet-class}} and the documentation of \code{eSet} (package Biobase) for the meaning of the other slots, which CountDataSet
   inherits from eSet (but which the present package does not use).
}
\value{
   an object of class CountDataSet
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
countsTable <- counts( makeExampleCountDataSet() )
cds <- newCountDataSet( countsTable, c( "A", "A", "A", "B", "B" ) )
}
