\name{subsetByGenes}
\alias{subsetByGenes}
\title{
    Making an ExonCountSet object from another one with a subset of its genes.
}
\description{
    Generates a smaller ExonCountSet object containing a subset of genes 
    from another ExonCountSet.
}
\usage{
subsetByGenes(ecs, genes)
}
\arguments{
  \item{ecs}{
      An ExonCountSet.
  }
  \item{genes}{
      Subset of geneIDs used to generate the subset ExonCountSet.
  }
}
\examples{
	data("pasillaExons", package="pasilla")
	ecs <- subsetByGenes(pasillaExons, sample(unique(geneIDs(pasillaExons)), 10))
}
