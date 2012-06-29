\name{ExonCountSet-class}
\docType{class}
\alias{ExonCountSet-class}

\title{"ExonCountSet", a container for exon count data }
\description{	This is the principal class of DEXSeq package. }
\section{Objects from the Class}{
Objects must be created with the function \code{\link{newExonCountSet}} (q.v.), alternatively
the user can call the function \code{\link{read.HTSeqCounts}}, which will call newExonCountset.
}
\section{Extends}{
Class \code{eSet} (package 'Biobase'), directly.
Class \code{VersionedBiobase} (package 'Biobase'), by class "eSet", distance 2.
Class \code{Versioned} (package 'Biobase'), by class "eSet", distance 3.
}
\note{ 
An ExonCountSet object stores the exon counts from high-throughput RNA sequencing experiments.
It is the principal object of the DEXSeq package. Some of the slots can be added by the user 
(see details in newExonCountSet documentation) or alternatively, the user can fill some of the slots
by using the HTSeq preprocessing steps and further calling \code{\link{read.HTSeqCounts}}, especially
those with the exon annotation data. The other slots will be filled with the analysis.

The ExonCountSet object contains a matrix of non-negative integers which represents sequence counts, 
with each column representing a sample and and each row a counting bin (i.e., an exon or part of an exon).  
In the phenoData, the object contains information about the samples,
e.g., size factors and design annotations are stored there. The user can also add more information about the
other properties of the samples.

An ExonCountSet object can be created just by providing a count matrix, and two vectors of gene and exon
identifiers of each of the rows in the matrix.  Nevertheless, the visualization plots included in DEXSeq
requires additional information about the exons (chromosome, strand, start, end).  This information can be added
directly after the creation of the ExonCountSet object.  If \code{\link{read.HTSeqCounts}} is called to create 
an ExonCountSet object, this information of the phenoData is inserted directly.

The columns for size factors (in phenoData), dispersion estimates, pvalue and padjust in the featureData
are filled later throughout the analysis, when the user calls \code{\link{estimateSizeFactors}}, \code{\link{estimateDispersions}}
\code{\link{fitDispersionFunction}}, and \code{\link{testForDEU}}.

}
\examples{
# See the vignette
}
