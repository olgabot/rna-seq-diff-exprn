setClass( "ExonCountSet", 
   contains = "eSet",
   representation = representation( 
      designColumns = "character",
      dispFitCoefs = "numeric",
      formulas = "list",
      annotationFile = "character"
   ),
   prototype = prototype( new( "VersionedBiobase",
      versions = c( classVersion("eSet"), ExonCountSet = "1.0.3" ) ) )
)


newExonCountSet <- function( countData, design, geneIDs, exonIDs, exonIntervals=NULL, transcripts=NULL){

   countData <- as.matrix( countData )
   if( any( round( countData ) != countData ) ){
      stop( "The countData is not integer." )}
   mode( countData ) <- "integer"

   if( is( design, "matrix" ) ){
      design <- as.data.frame( design )}

#   if(!(is(design, "data.frame") || is(design, "AnnotatedDataFrame"))){
 #     countData <- countData[,order(as.character(design))]
  #    design <- design[order(as.character(design))]
 #  }else{
 #     countData <- countData[,order(design$condition)]
 #     design <- design[order(design$condition),]
 #  }

   phenoData <- annotatedDataFrameFrom( countData, byrow=FALSE )
   featureData <- annotatedDataFrameFrom( countData, byrow=TRUE )
      
   phenoData$sizeFactor <- rep( NA_real_, ncol(countData) )
   varMetadata( phenoData )[ "sizeFactor", "labelDescription" ] <- "size factor (relative estimate of sequencing depth)"

   geneIDs <- as.factor( geneIDs )
   if( length(geneIDs) != nrow(countData) )
      stop( "geneIDs must be of the same length as the number of columns in countData")

   featureData$geneID <- geneIDs
   varMetadata( featureData )[ "geneID", "labelDescription" ] <- "ID of gene to which the exon belongs"

   exonIDs <- as.character( exonIDs )
   if( length(exonIDs) != nrow(countData) )
      stop( "exonIDs must be of the same length as the number of columns in countData")

   featureData$exonID <- exonIDs
   varMetadata( featureData )[ "exonID", "labelDescription" ] <- "exon ID (unique only within a gene)"

   if( is.null(exonIntervals) ){
      exonIntervals <- data.frame(
         chr    = rep( NA_character_, nrow( countData ) ), 
         start  = rep( NA_integer_,   nrow( countData ) ), 
         end    = rep( NA_integer_,   nrow( countData ) ), 
         strand = rep( NA_character_, nrow( countData ) ) ) }

   featureData$testable <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "testable", "labelDescription" ] <- "slot indicating if an exon should be considered in the test"

   featureData$dispBeforeSharing <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "dispBeforeSharing", "labelDescription" ] <- "exon dispersion (Cox-Reid estimate)"

   featureData$dispFitted <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "dispFitted", "labelDescription" ] <- "Fitted mean-variance estimate alpha(mu)=alpha1/mu + alpha0"

   featureData$dispersion <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "dispersion", "labelDescription" ] <- "maximum value between the exon dispersion before sharing and the fitted dispersion estimate"

   featureData$pvalue <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "pvalue", "labelDescription" ] <- "p-value from testForDEU"

   featureData$padjust <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "padjust", "labelDescription" ] <- "BH adjusted p-value"

   exonIntervals <- as.data.frame( exonIntervals )
   
   # in case it was a GRanges object before, change the colname:
   if( "seqnames" %in% colnames(exonIntervals) ){
      colnames(exonIntervals)[ colnames(exonIntervals) == "seqnames" ] <- "chr"   }
   
   if( !all( c( "chr", "start", "end", "strand" ) %in% colnames(exonIntervals) ) ){
      stop( "exonIntervals must be a data frame with columns 'chr', 'start', 'end', and 'strand'." )}

   if(is.null(transcripts)){
      transcripts <- rep(NA_character_, nrow( countData ) )}

   if(!is(transcripts, "character")){
      stop("transcript information must be a character vector")}
   
   featureData$chr    <- factor( exonIntervals$chr )
   featureData$start  <- exonIntervals$start
   featureData$end    <- exonIntervals$end
   featureData$strand <- factor( exonIntervals$strand )
   featureData$transcripts <- transcripts
   varMetadata( featureData )[ "chr",    "labelDescription" ] <- "chromosome of exon"
   varMetadata( featureData )[ "start",  "labelDescription" ] <- "start of exon"
   varMetadata( featureData )[ "end",    "labelDescription" ] <- "end of exon"
   varMetadata( featureData )[ "strand", "labelDescription" ] <- "strand of exon"
   varMetadata( featureData )[ "transcripts", "labelDescription" ] <- "transcripts in which this exon is contained"
   
   if( is( design, "data.frame" ) || is( design, "AnnotatedDataFrame" ) ) {
      stopifnot( nrow( design ) == ncol( countData ) )
      design <- as( design, "AnnotatedDataFrame" )
      dimLabels(design) <- dimLabels(phenoData)
      rownames( pData(design) ) <- rownames( pData(phenoData) )
      phenoData <- combine( phenoData, design )
      rvft <- c( `_all` = NA_character_ )
      designColumns <- varLabels(design)
   } else {
      design <- factor( design, levels=unique(design))
      stopifnot( length( design ) == ncol( countData ) )
      phenoData$`condition` <- factor( design )
      varMetadata( phenoData )[ "condition", "labelDescription" ] <- "experimental condition, treatment or phenotype"
      designColumns <- "condition"
   }
   ecs <- new( "ExonCountSet",
      assayData = assayDataNew( "environment", counts=countData ),
      phenoData = phenoData, 
      featureData = featureData,
      designColumns = designColumns,
      dispFitCoefs = c( NA_real_, NA_real_ )
      )
   ecs
}

setValidity( "ExonCountSet", function( object ) {

   if( !all( object@designColumns %in% names(pData(object)) ) )
      return( "Not all designColumns appear in phenoData." )

   if( ! "sizeFactor" %in% names(pData(object)) )
      return( "phenoData does not contain a 'sizeFactor' column.")
   if( ! is( pData(object)$`sizeFactor`, "numeric" ) )
      return( "The 'sizeFactor' column in phenoData is not numeric." )

   if( ! "geneID" %in% names(fData(object)) )
      return( "featureData does not contain a 'geneID' column.")
   if( ! is( fData(object)$geneID, "factor" ) )
      return( "The 'geneID' column in fData is not a factor." )

   if( ! "exonID" %in% names(fData(object)) )
      return( "featureData does not contain an 'exonID' column.")
   if( ! is( fData(object)$exonID, "character" ) )
      return( "The 'exonID' column in fData is not a character vector." )

   if( ! "chr"  %in% names(fData(object)) )
      return( "featureData does not contain a 'chr' column.")
   if( ! is( fData(object)$chr, "factor" ) )
      return( "The 'chr' column in fData is not a factor." )

   if( ! "start"  %in% names(fData(object)) )
      return( "featureData does not contain a 'start' column.")
   if( ! is( fData(object)$start, "integer" ) )
      return( "The 'start' column in fData is not integer." )

   if( ! "end"  %in% names(fData(object)) )
      return( "featureData does not contain a 'end' column.")
   if( ! is( fData(object)$end, "integer" ) )
      return( "The 'end' column in fData is not integer." )

   if( ! "strand"  %in% names(fData(object)) )
      return( "featureData does not contain a 'strand' column.")
   if( ! is( fData(object)$strand, "factor" ) )
      return( "The 'strand' column in fData is not a factor." )
#   if( any( sort( levels(fData(object)$strand)) != c( "-", "+" ) ) )
#      return( "The 'strand' column in fData does not have the levels '+' and '-'." )   ### strange reason, make pasilla check crash only in the check, not sourcing exactly the same code
   if( !is(fData(object)$dispersion, "numeric")){
      return( "The 'dispersion' is not numeric")}
   if( !is(fData(object)$dispFitted, "numeric")){
      return( "The 'dispFitted' is not numeric")}
   if( !is(fData(object)$dispBeforeSharing, "numeric")){
      return( "The 'dispBeforeSharing' column is not numeric")}
   if( !is(fData(object)$pvalue, "numeric")){
      return( "The 'pvalue' values are not numeric")}
   if( !is(fData(object)$padjust, "numeric")){
      return( "The 'padjust' values are not numeric")}
   if( !is.integer( assayData(object)[["counts"]] ) )
      return( "The count data is not in integer mode." )

   if( any( assayData(object)[["counts"]] < 0 ) )
      return( "The count data contains negative values." )

   if( length( object@dispFitCoefs ) != 2 )
      return( "dispFitCoefs is not a vector of length 2." )

   TRUE
} )

setMethod("counts", signature(object="ExonCountSet"),
  function( object, normalized=FALSE) {
    cds <- object
    if(!normalized){
      assayData(cds)[["counts"]]
    } else {
      if(any(is.na( sizeFactors(cds)))) {
         stop( "Please first calculate size factors or set normalized=FALSE")
      } else {
         t(t( assayData( cds )[["counts"]] ) / sizeFactors(cds) )
      }
   }
})

setReplaceMethod("counts", signature(object="ExonCountSet", value="matrix"),
  function( object, value ) {
   cds <- object
   assayData(cds)[[ "counts" ]] <- value
   validObject(cds)
   cds
})   

setMethod("sizeFactors",  signature(object="ExonCountSet"),
  function(object) {
   cds <- object
   sf <- pData(cds)$sizeFactor
   names( sf ) <- colnames( counts(cds) )
   sf
})

setReplaceMethod("sizeFactors",  signature(object="ExonCountSet", value="numeric"),
  function(object, value ) {
   cds <- object
   pData(cds)$sizeFactor <- value
   validObject( cds )
   cds
})


setMethod("design", signature(object="ExonCountSet"),
   function( object, drop=TRUE, asAnnotatedDataFrame=FALSE ) {
      cds <- object
      if( asAnnotatedDataFrame )
         return( phenoData(cds)[, cds@designColumns ] )
      ans <- pData(cds)[, cds@designColumns, drop=FALSE ]
      if( ncol(ans) == 1 && drop ) {
         ans <- ans[,1]
         names(ans) <- colnames( counts(cds) ) }
      else
         rownames( ans ) <- colnames( counts(cds) )
      ans
})

setReplaceMethod("design", signature(object="ExonCountSet"),
      function( object, value ) {
      cds <- object
      ## Is it multivariate or just a vector?
      if( ncol(cbind(value)) > 1 )
         value <- as( value, "AnnotatedDataFrame" )
      else {
         value <- new( "AnnotatedDataFrame",
            data = data.frame( condition = value ) )
         varMetadata( value )[ "condition", "labelDescription" ] <-
            "experimental condition, treatment or phenotype" }

      rownames( pData(value) ) <- rownames( pData(cds) )
      dimLabels( value ) <- dimLabels( phenoData(cds) )
      phenoData(cds) <- combine( 
         phenoData(cds)[ , !( colnames(pData(cds)) %in% cds@designColumns ), drop=FALSE ], 
         value )
      cds@designColumns <- colnames( pData(value) )   
      validObject(cds)
      cds
})

#setMethod("design",       signature(cds="ExonCountSet"), .design)
#setMethod("design<-",     signature(cds="ExonCountSet"), `.design<-`)
#setMethod("conditions",   signature(cds="ExonCountSet"), .design)
#setMethod("conditions<-", signature(cds="ExonCountSet"), `.design<-`)

geneIDs <- function( ecs ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   g <- fData(ecs)$geneID
   names(g) <- rownames( counts(ecs) )
   g
}

`geneIDs<-` <- function( ecs, value ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   fData(ecs)$geneID <- value
   validObject(ecs)
   ecs
}
      
exonIDs <- function( ecs ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   g <- fData(ecs)$exonID
   names(g) <- rownames( counts(ecs) )
   g
}

`exonIDs<-` <- function( ecs, value ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   fData(ecs)$exonID <- value
   validObject(ecs)
   ecs
}
            

subsetByGenes <- function( ecs, genes ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   stopifnot( all( genes %in% levels(geneIDs(ecs)) ) )
   ecs2 <- ecs[ as.character(geneIDs(ecs)) %in% genes, ]
   ecs2
}
   
geneCountTable <- function( ecs ) {
   stopifnot( is( ecs, "ExonCountSet" ) )
   do.call( rbind, 
      tapply( 1:nrow(ecs), geneIDs(ecs), function(rows) 
         colSums( counts(ecs)[rows,,drop=FALSE] ) ) )
}

DEUresultTable <- function(ecs)
{
   result <- data.frame(
      geneID=geneIDs(ecs), 
      exonID=exonIDs(ecs), 
      dispersion=featureData(ecs)$dispersion, 
      pvalue=fData(ecs)$pvalue, 
      padjust=fData(ecs)$padjust,
      meanBase=rowMeans(counts(ecs, normalized=TRUE)))

   extracol <- regexpr("log2fold", colnames(fData(ecs)))==1
   if(any(extracol)){
      w <- which(extracol)
      result <- data.frame(result, fData(ecs)[,w])
      colnames(result)[7:(6+length(w))] <- colnames(fData(ecs))[w]
   }
   result
}
