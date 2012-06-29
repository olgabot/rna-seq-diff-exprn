\name{CountDataSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CountDataSet-class}

\title{Class "CountDataSet" -- a container for count data from HTS experiments }
\description{	This is the main class for the present package. }
\section{Objects from the Class}{
Objects should be created with calls to \code{\link{newCountDataSet}} (q.v.).
}
\section{Extends}{
Class \code{eSet} (package 'Biobase'), directly.
Class \code{VersionedBiobase} (package 'Biobase'), by class "eSet", distance 2.
Class \code{Versioned} (package 'Biobase'), by class "eSet", distance 3.
}
\note{ 
Note: This is a summary for reference. For an explanation of the actual usage, see the vignette.

A CountDataSet object stores counts from an HTS data set and offers
further slots which are populated during the analysis.

After creation with \code{\link{newCountDataSet}}, a CountDataSet typically contains 
a count table, i.e., a matrix of integer data, that is accessible with the accessor function \code{\link{counts}}.
Each row of the matrix corresponds to a gene (or binding region, or the like), and
each colum to an experimental sample. The experimental conditions of the samples
are stored in a factor (with one element for each row of the counts matrix), which
can be read with the accessor function \code{\link{conditions}}.

In the following analysis steps, further data slots are populated. First, the size factors
can be estimated with \code{\link{estimateSizeFactors}}, which are afterwards accessible 
via \code{\link{sizeFactors}}. Then, the dispersions (variance fits) are estimated with \code{\link{estimateDispersions}}. The resulting estimates are stored in phenoData columns, accessible
via \code{pData}, with the column names staring with \code{disp_}. The intermediate steps
of the fit are stored in the environment-values slot \code{fitInfo} (see \code{\link{estimateDispersions}} 
for details).

Internally, the mentioned data is stored in slots as follows: 

As \code{CountDataSet} is derived from \code{eSet}, it has a \code{phenoData} slot
which allows to store sample annotation. This is used to store the factor with the conditions,
as a data frame column named \code{condition}, and to store the size factors, as an numeric data 
frame column named \code{sizeFactor}. If the user creates an object with multivariate design,
i.e., passes a data frame instead of a factor for \code{conditions}, this data frame's
columns are placed in the \code{phenoData} slot instead of the \code{condition} column.
Furthermore, the function \code{\link{estimateDispersions}} adds columns with the dispersion
values to be used by \code{\link{nbinomTest}} and \code{\link{fitNbinomGLMs}}. These
columns have names starting with \code{disp_}.

The user may add further columns to the \code{phenoData} AnnotatedDataFrame.

The counts table is stored in the \code{eSet}'s \code{assayData} locked environment with 
the name \code{counts}.

The slot \code{dispInfo} is an environment containing lists, one for each set of estimated
dispersion values and the slot \code{dispTable} (with accessor \code{\link{dispTable}}
shows the assignment of conditions to dispersion estimates. See \code{\link{estimateDispersions}}
}

\examples{
# See the vignette
}

