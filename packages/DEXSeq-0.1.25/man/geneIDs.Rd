\name{geneIDs}
\alias{geneIDs}
\alias{geneIDs<-}
\title{
   Accessor for the geneIDs in an ExonCountSet object.
}
\description{
   This function is an accessor for the gene identifiers for each of the rows in the count table.
}
\usage{
geneIDs(ecs)
geneIDs(ecs) <- value
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
  \item{value}{
     An factor of gene identifiers, one for each of the rows of the count data.
   }
}
\examples{
data("pasillaExons", package="pasilla")
head( geneIDs( pasillaExons ) )
}   
