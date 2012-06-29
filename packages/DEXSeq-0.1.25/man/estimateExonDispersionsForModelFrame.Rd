\name{estimateExonDispersionsForModelFrame}
   \alias{estimateExonDispersionsForModelFrame}
   \title{Estimates exon dispersions}
   \description{
       This function calculates the individual dispersions for each
       counting bins for a single gene. It takes as input  a 
       model frame generated from the function \code{\link{modelFrameForGene}}.
	}
\usage{
estimateExonDispersionsForModelFrame(modelFrame, formula=NULL, mm=NULL, 
                                                 muhat=NULL, initialGuess=0.01)
}
\arguments{
     \item{modelFrame}{
        Model frame provided by the function \code{\link{modelFrameForGene}}.
     }
     \item{formula}{
        Formula for the glm used to estimate the dispersions. The factors in the formula
        must be present in the column names of the model frame. If it is left NULL, the default 
        formula used is "count ~ sample + condition * exon".
     }
     \item{mm}{
        A model matrix for the model frame. If NULL, a model matrix will 
        be created from the parameters "formula" and "modelFrame".
     }
     \item{muhat}{
        Initial values for the coefficients in the optimization.  
        If NULL, initial values will be calculated using with the dispersion
        value given by the parameter "initialGuess".
     }
     \item{initialGuess}{
        An initial guess of the dispersions to initiate the optimization.
     }
}
\value{
	A vector of exon dispersions.
}
\examples{
	data("pasillaExons", package="pasilla")
	pasillaExons <- estimateSizeFactors( pasillaExons )
	estimateExonDispersionsForModelFrame(modelFrameForGene(pasillaExons, "FBgn0085442"))
}
\seealso{
        \code{\link{estimateDispersions}}
}
