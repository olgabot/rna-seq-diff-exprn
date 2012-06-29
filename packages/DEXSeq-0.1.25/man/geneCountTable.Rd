\name{geneCountTable}
\alias{geneCountTable}
\title{
   Makes a count table for genes.
}
\description{
   This function returns a count table  where each row is a gene and each column is a sample, by adding up
   the values for each gene's individual counting bins.
}
\usage{
geneCountTable(ecs)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
}
}
\seealso{
   \code{DESeq}
}   
\examples{
data("pasillaExons", package="pasilla")
head(geneCountTable(pasillaExons))
}   

