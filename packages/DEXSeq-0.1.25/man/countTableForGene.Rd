\name{countTableForGene}
\alias{countTableForGene}
\title{
    Count table for a given geneID.
}
\description{
   This function returns a matrix of non negative integers containing 
   a count table for a specified geneID from an ExonCountSet object.  
   The count table contains one row for every counting bin of the gene and 
   a column for every sample.
}
\usage{
countTableForGene(ecs, geneID, normalized=FALSE, withDispersion=FALSE)
}
\arguments{
\item{ecs}{
   An ExonCountSet.
}
\item{geneID}{
   A geneID to get the count table.
}
\item{normalized}{
   If TRUE, the raw counts will be normalized by the size factors.
}
\item{withDispersion}{
   If TRUE, an extra column with the dispersion estimate used in the 
   test will added to the count table.
}
}

\seealso{
\code{\link{estimateSizeFactors}}
}
\examples{
   data("pasillaExons", package="pasilla")
   pasillaExons <- estimateSizeFactors( pasillaExons )
   countTableForGene(pasillaExons, "FBgn0085442", normalized=FALSE)
}
