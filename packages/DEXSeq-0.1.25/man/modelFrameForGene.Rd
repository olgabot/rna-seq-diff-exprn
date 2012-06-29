\name{modelFrameForGene}
\alias{modelFrameForGene}
\title{
   Makes the model frame for a geneID.
}
\description{
   Creates a data frame containing the model frame for a gene with the columns 
   sample, exon, size factors, their respective counts and the design annotation.
}
\usage{
modelFrameForGene(ecs, geneID, onlyTestable=FALSE)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
  \item{geneID}{
     A gene identificator contained in the ExonCountSet object.
  }
  \item{onlyTestable}{
     If TRUE, only the testable exons will be included in the model frame.  
     Check \code{fData$testable} for more information.
  }
}
\value{
     A data frame containing the model frame for a gene.
}
\examples{
data("pasillaExons", package="pasilla")
modelFrameForGene(pasillaExons, "FBgn0085442")
}

