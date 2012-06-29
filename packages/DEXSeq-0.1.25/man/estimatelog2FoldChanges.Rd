\name{estimatelog2FoldChanges}
   \alias{estimatelog2FoldChanges}
   \title{Fold changes (log2) from the fitted expression values in the GLM.}
   \description{
      This function calculates the fold changes (on log2 scale) between the different
      conditions. It calculates them from the coefficients of a GLM
      that fits the read counts to a variable of the experimental design specified 
      by the user (see below, parameter "fitExpToVar").
   }
   \usage{
estimatelog2FoldChanges(ecs, fitExpToVar="condition", 
          nCores=1, quiet=FALSE, file="")
   }
   \arguments{
      \item{ecs}{
         An ExonCountSet object.
      }
      \item{fitExpToVar}{
         A variable contained in \code{design(ecs)}. The expression 
         values will be fitted to this variable using the the formula 
         "count ~ sample + fitExpToVar * exon".
      }
      \item{nCores}{
         Number of CPU cores to be used to estimate the dispersions.  
         The \code{multicore} package need to be loaded beforehands to
         parallelize over several cores.
      }
      \item{quiet}{
         If TRUE, no progress report is shown.  In case the session is 
         not an interactive and progress report is wanted, add a file name below.
      }
      \item{file}{
         A file name to write the progress reports. If file="", output
         will be written to the standard output connection.
      }
  }
  \examples{
     \dontrun{
         data("pasillaExons", package="pasilla")
         pasillaExons <- estimateSizeFactors( pasillaExons )
         pasillaExons <- estimateDispersions( pasillaExons )
         pasillaExons <- fitDispersionFunction( pasillaExons ) 
         pasillaExons <- estimatelog2FoldChanges( pasillaExons )
     }
  }
\keyword{ExonCountSet}
