\name{nbinomTestForMatrices}
\Rdversion{1.1}
\alias{nbinomTestForMatrices}
\title{
   Perform row-wise tests for differences between the base means of two count matrices.
}
\description{
   This function is called by \code{\link{nbinomTest}}. Call it directly only
   if the S4 interface is unsuitable for your task.
}
\usage{
nbinomTestForMatrices(countsA, countsB, sizeFactorsA, sizeFactorsB, dispsA, dispsB )
}
\arguments{
  \item{countsA}{ 
     A matrix of counts, where each column is a replicate
}
  \item{countsB}{
     Another matrix of counts, where each column is a replicate}
  \item{sizeFactorsA}{
     Size factors for the columns of the matrix 'countsA'
}
  \item{sizeFactorsB}{
     Size factors for the columns of the matrix 'countsB'
}
  \item{dispsA}{
     The dispersions for 'countsA', a vector with one value per gene
}
  \item{dispsB}{
     The same for 'countsB'
}
}
\details{
   See the vignette for an exact description of the null hypothesis tested.
}
\value{
   A vector of unadjusted p values, one for each row in the counts matrices.
}

\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds, method="per-condition" )
colsA <- conditions(cds) == "A"
colsB <- conditions(cds) == "B"
bmvA <- getBaseMeansAndVariances( counts(cds)[,colsA], sizeFactors(cds)[colsA] )
bmvB <- getBaseMeansAndVariances( counts(cds)[,colsB], sizeFactors(cds)[colsB] )
pvals <- nbinomTestForMatrices( 
   counts(cds)[,colsA],
   counts(cds)[,colsB],
   sizeFactors(cds)[colsA],
   sizeFactors(cds)[colsB],
   fitInfo(cds,"A")$dispFunc( rowMeans( counts( cds, normalized=TRUE ) ) ),
   fitInfo(cds,"B")$dispFunc( rowMeans( counts( cds, normalized=TRUE ) ) ) )
names( pvals ) <- row.names( counts(cds) )
head( pvals )

# This here should give the same results:
head( nbinomTest( cds, "A", "B" )$pval )
}
