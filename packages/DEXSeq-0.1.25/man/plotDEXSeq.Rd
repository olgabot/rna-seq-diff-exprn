\name{plotDEXSeq}
     \alias{plotDEXSeq}
     \title{Visualization of the fitted expression, fitted splicing or the normalized counts.}
     \description{
	The function provides a plot to visualize read count data, the fitted 
        expression, fitted splicing and the results of the test in 
        \code{\link{testForDEU}}.  The fitted values are obtained from fitting
        the counts values to a certain condition from the design annotation of
        the glm. See fitExpToVar parameter.
     }
     \usage{
plotDEXSeq(ecs, geneID, FDR=0.1, fitExpToVar="condition", 
           norCounts=FALSE, expression=TRUE, splicing=FALSE, 
           displayTranscripts=FALSE, names=FALSE, legend=FALSE, 
           color=NULL, color.samples=NULL, ...)
     }
     \arguments{
       \item{ecs}{
          An ExonCountSet object.
       }
       \item{geneID}{
          ID of the gene to visualize.
       }
       \item{FDR}{
          A false discovery rate used to indicate the significant exons.
       }
       \item{fitExpToVar}{
          A variable contained in the design annotation of the \code{ExonCountSet}, 
          the expression values will be fitted to this variable using the formula
          \code{count ~ fitExpToVar * exon} using a model frame obtained from the function
          \code{\link{modelFrameForGene}}.  
       }
       \item{norCounts}{
          If TRUE, provides a plot of the counts normalized by the size factors.
       }
       \item{expression}{
          If TRUE, the function plots the fitted EXPRESSION estimates from the glm 
          regression.
       }
       \item{splicing}{
          If TRUE, the function plots the fitted SPLICING estimates from the glm 
          regression.
       }
       \item{displayTranscripts}{
          If TRUE, the transcripts are displayed in the plot.
       }
       \item{names}{
          If TRUE, the names of the transcripts are shown.
       }
       \item{legend}{
          If TRUE, a legend is displayed.
       }
       \item{color}{
          A vector of colors for each of the levels of the factor in the design of 
          the ExonCountSet object indicated by "fitExpToVar".
       }
       \item{color.samples}{
          A vector of colors for each of the samples. If NULL, the colors of each 
          sample will be assigned according to its corresponding level from "fitExpToVar". 
          This option is useful to visualize complex experimental designs.
       }
       \item{...}{
          Further graphical parameters (see \code{par}).
       }
     }
     \seealso{
       graphics, \code{segments}
     }
     \examples{
	\dontrun{
           data("pasillaExons", package="pasilla")
	   pasillaExons <- estimateSizeFactors(pasillaExons)
	   pasillaExons <- estimateDispersions(pasillaExons)
           pasillaExons <- fitDispersionFunction( pasillaExons )
	   plotDEXSeq(pasillaExons, "FBgn0085442")
        }
     }
\keyword{ExonCountSet}
