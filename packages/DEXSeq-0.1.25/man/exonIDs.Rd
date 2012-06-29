\name{exonIDs}
\alias{exonIDs}
\alias{exonIDs<-}
\title{
   Accessor for the exonIDs in an ExonCountSet object.
}
\description{
   This function is an accessor for the exon identifiers for each of the rows in the count table.
   Note that each exon ID identifies, strictly speaking, not an exon but a counting bin,
   which may well be just part of an exon.
   Make sure that the exon IDs are ordered alphanumerically in the gene.
}
\usage{
exonIDs(ecs)
exonIDs(ecs) <- value
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
   }
  \item{value}{
     A vector of exon counting bin identifiers, one for each of the rows of the count data.
   }
}
\examples{
	library(DEXSeq)
   	data("pasillaExons", package="pasilla")
	exonIDs(pasillaExons)
}
