plotDEXSeq <- function(ecs, geneID, FDR=0.1, fitExpToVar="condition", norCounts=FALSE, expression=TRUE, splicing=FALSE, displayTranscripts=FALSE, names=FALSE, legend=FALSE, color=NULL, color.samples=NULL, ...)
{
   stopifnot(is(ecs, "ExonCountSet"))
   if(any(is.na(sizeFactors(ecs)))){
      stop("Please estimate sizeFactors first\n")}
   if(!fitExpToVar %in% ecs@designColumns){
      stop("fitExpToVar parameter is not in the design columns, double check ecs@designColumns")}
   if(sum(is.na(featureData(ecs)$dispersion))==nrow(counts(ecs))){
      stop("No dispersion parameters found, first call function fitDispersionFunction...\n")}
   op <- sum(c(expression, splicing, norCounts))
   if(op == 0){
      stop("Please indicate what would you like to plot\n")}

   rt<-which(featureData(ecs)$geneID==geneID)
   count <- countTableForGene(ecs, geneID, normalized=TRUE)
   if(sum(count) == 0){
      warning("No read counts falling in this gene, there is nothing to plot.")
      return()}
   if(FDR>1|FDR<0){
      stop("FDR has to be a numeric value between 0 - 1")}

   rango <- 1:length(rt)
   intervals<-(0:nrow(count))/nrow(count)
   numcond<-length(unique(design(ecs, drop=FALSE)[[fitExpToVar]]))
   numexons<-nrow(count)
   each <- featureData(ecs)$padjust[rt]
   #exoncol<-ifelse(each<=FDR, "#8B0000", "dark green")
   #exoncol[is.na(exoncol)]<-"black"
   #colorlines <- ifelse(each<=FDR, "#FF000060", "lightgrey")
   exoncol<-ifelse(each<=FDR, "#F219ED", "#CCCCCC")
   exoncol[is.na(exoncol)]<-"white"
   colorlines <- ifelse(each<=FDR, "#F219ED60", "#B3B3B360")   # vertical dashed lines
   colorlines[is.na(colorlines)] <- "#B3B3B360"
   colorlinesB <- ifelse(each<=FDR, "#9E109B", "#666666")  # slanted solid lines
   colorlinesB[is.na(colorlinesB)] <- "#666666"


   ################## DETERMINE THE LAYOUT OF THE PLOT DEPENDING OF THE OPTIONS THE USER PROVIDES ###########
   if(!any(is.na(featureData(ecs)$start))){
      sub<-data.frame(start=featureData(ecs)$start[rt], end=featureData(ecs)$end[rt], chr=featureData(ecs)$chr[rt], strand=featureData(ecs)$strand[rt])
      rel<-(data.frame(sub$start, sub$end))-min(sub$start)	
      rel<-rel/max(rel[,2])
      if(displayTranscripts==TRUE & !is.null(featureData(ecs)$transcripts)){
         transcripts <- sapply(featureData(ecs)$transcripts[rt], function(x){strsplit(x, ";")})
         trans <- Reduce(union, transcripts)
         if(length(trans) > 40){
            warning("This gene contains more than 40 transcripts annotated, only the first 40 will be plotted\n")
         }
         mat <- 1:(3+min(length(trans), 40)) ## max support from transcripts is 45, which seems to be the max for the layout supported by graphics
         hei<-c(8, 1, 1.5, rep(1.5, min(length(trans), 40)))
      }else{
         mat<-1:3
         hei<-c(5, 1, 1.5)
      }
      if(op > 1){
         hei <- c(rep(hei[1], op-1), hei)
         mat <- c(mat, (length(mat)+(1:length(op))))
      }
      hei <- c(hei, .2)
      mat <- c(mat, length(mat)+1)
      layout(matrix(mat), heights=hei)
      par(mar=c(2, 4, 4, 2))
   }else if(op > 1){
		par(mfrow=c(op,1))
   }
   
   ####### DETERMINE COLORS, IF THE USER DOES NOT PROVIDE ONE PER SAMPLE THE COUNT WILL OBTAIN THEM CORRESPONDING TO THEIR DESIGN ####
   ##### determine colors if not provided by user ######
   if(is.null(color)){
      color<-rgb(colorRamp(c("#D7191C", "#FFFFBF", "#2B83BA"))(seq(0, 1, length.out=numcond)), maxColorValue=255, alpha=175)
   }
   names(color) <- sort(levels(design(ecs, drop=FALSE)[[fitExpToVar]]))

   
   if(expression){
      es <- fitAndArrangeCoefs( ecs, geneID, frm=as.formula(paste("count ~", fitExpToVar,  "* exon")) )
      if(is.null(es)){
          warning(sprintf("glm fit failed for gene %s", geneID))
          return()
      }
      coeff <- as.matrix( t( getEffectsForPlotting(es, averageOutExpression=FALSE, groupingVar=fitExpToVar) ) )
      coeff <- exp(coeff)
      ylimn <- c(0, max(coeff, na.rm=TRUE))
      coeff <- vst( coeff, ecs )
      drawPlot(matr=coeff, ylimn, ecs, intervals, rango, textAxis="Fitted expression", rt=rt, color=rep(color[sort(levels(design(ecs, drop=FALSE)[[fitExpToVar]]))], each=numexons), colorlines=colorlines, ...)
   }

   if(splicing){
      coeff <- as.matrix( t( getEffectsForPlotting( fitAndArrangeCoefs( ecs, geneID, frm=as.formula(paste("count ~", fitExpToVar,  "* exon")) ), averageOutExpression=TRUE, groupingVar=fitExpToVar) ) )
      coeff <- exp(coeff)
      ylimn <- c(0, max(coeff, na.rm=TRUE))
      coeff <- vst( coeff, ecs )
      drawPlot(matr=coeff, ylimn, ecs, intervals, rango, textAxis="Fitted splicing", rt=rt, color=rep(color[sort(levels(design(ecs, drop=FALSE)[[fitExpToVar]]))], each=numexons), colorlines=colorlines, ...)
   }

   if(norCounts){
      ylimn <- c(0, max(count, na.rm=TRUE))
      count <- vst( count, ecs )
      if(is.null(color.samples)){
         colorcounts <- rep(color[as.character(design(ecs, drop=FALSE)[[fitExpToVar]])], each=numexons)
      }else{
         colorcounts <- rep(color.samples, each=numexons)
      }
         drawPlot(matr=count, ylimn, ecs, intervals, rango, textAxis="Normalized counts", rt=rt, color=colorcounts, colorlines=colorlines, ...)
   }

	########### plot the gene model ########## just if transcript information available
   if( !any(is.na(featureData(ecs)$start)) ){
      par(mar=c(0, 4, 0, 2))
      plot.new()
      segments(apply((rbind(rel[rango,2], rel[rango, 1])), 2, median), 0, apply(rbind(intervals[rango], intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2)), 2, median), 1, col=colorlinesB)
      par(mar=c(1.5, 4, 0, 2))
      drawGene(min(sub$start), max(sub$end), tr=sub, rango, exoncol=exoncol, names, trName="Gene model", cex=0.8)
      if(!is.null(featureData(ecs)$transcripts)){
      ##### plot the transcripts #######
         if(displayTranscripts){
            for(i in 1:min(length(trans), 40)){
               logicexons <- sapply(transcripts, function(x){length(which(x==trans[i]))})
               tr<-data.frame(featureData(ecs)$start[rt][logicexons==1], featureData(ecs)$end[rt][logicexons==1])
               drawGene(min(sub$start), max(sub$end), tr=tr, rango, exoncol=NULL, names, trName=trans[i], cex=0.8)
            }
         }
      }
      axis(1, at=round(seq(min(sub$start), max(sub$end), length.out=10)), labels=round(seq(min(sub$start), max(sub$end), length.out=10)), pos=0, lwd.ticks=0.2, padj=-0.7, ...)   ########## genome axis
   }
   if(legend){
      mtext(paste(geneID, unique(featureData(ecs)$strand[rt])), side=3, adj=0.25, padj=1.5, line=0, outer=TRUE, cex=1.5)
      posforlegend <- seq(.7, .9, length.out=numcond)
      for(i in 1:length(color)){
         mtext(names(color[i]), side=3, adj=posforlegend[i], padj=1.5, line=0, outer=TRUE, col=color[i], ...)
      }
   }else{
      mtext(paste(geneID, unique(featureData(ecs)$strand[rt])), side=3, adj=0.5, padj=1.5, line=0, outer=TRUE, cex=1.5)
   }
}

###################
#VST FUNCTION
###################
vst <- function(x, ecs)
{
( 2 / ( sqrt(ecs@dispFitCoefs[1]) ) ) * 
   log( 2 * ecs@dispFitCoefs[1] * sqrt(x) + 
        2 * sqrt( ecs@dispFitCoefs[1] * ( ecs@dispFitCoefs[2] + 1 + ecs@dispFitCoefs[1] * x ) ) ) -
( 2 / ( sqrt(ecs@dispFitCoefs[1]) ) ) * 
   log( 2 * sqrt( ecs@dispFitCoefs[1] * ( ecs@dispFitCoefs[2] + 1 ) ) )
}

####################
#FUNCTION TO MAKE THE AXIS OF THE VST VALUES
####################

makevstaxis <- function(min, ylimn, ecs, ...)
{
   minlog10 <- floor( log10( 1/ncol(counts(ecs)) ) )
   maxlog10 <- ceiling( log10( ylimn[2] ) )
   ticks <- 10^seq( minlog10, maxlog10 )
   decade_lengths <- ( vst(ticks, ecs)[ 2 : length(ticks) ] - vst(ticks, ecs)[ 1 : (length(ticks)-1) ] ) /
      ( vst( ylimn[2], ecs) - vst( ylimn[1], ecs) )
#   ticks <- c( 0, ticks[ min( which( decade_lengths > .1 ) ) : length(ticks) ] )
   fromw <- pmax( which( decade_lengths > .1 ), 0 )[1]
   ticks <- ticks[ fromw : length(ticks) ]
   axis( 2, at=vst(c(0, ticks), ecs), labels=c("",ticks), las=2, pos=0, ...)

   for( i in minlog10 : (maxlog10-1) ) {
      decade_length <- ( vst( 10^(i+1), ecs) - vst( 10^i, ecs) ) / ( vst( ylimn[2], ecs) - vst( ylimn[1], ecs) )
      if( decade_length > .1 ) {
         axis( 2, at = vst( 1:9 * 10^i, ecs), labels = FALSE, las=2, tcl=-.25, pos=0, ...)
      }
      if( decade_length > .4 & decade_length <= .6) {
         axis( 2, at = vst( c(2,3,5) * 10^i, ecs ), labels = ( c(2,3,5) * 10^i ), las=2, tcl=-.25, pos=0, ...)
      } else if( decade_length > .6 ) {
         axis( 2, at = vst( c(1.5,2:9) * 10^i, ecs ), labels = ( c(1.5,2:9) * 10^i ), las=2, tcl=-.25, pos=0, ...)
      }
   }
#   axis( 2, at=vst(0, ecs), labels=FALSE, las=2, pos=0, ...)
}



#######################
#FUNCTION TO WRITE THE PLOTS:
#######################
drawPlot <- function(matr, ylimn, ecs, intervals, rango, fitExpToVar, numexons, textAxis, rt, color, colorlines, ...)
{
   plot.new()
   plot.window(xlim=c(0, 1), ylim=c(0, max(matr)))
   makevstaxis(1/ncol(matr), ylimn, ecs, ...)
   intervals<-(0:nrow(matr))/nrow(matr)
   middle <- apply(cbind(intervals[rango], (intervals[rango+1]-((intervals[rango+1])-intervals[rango])*0.2)), 1, median)
   matr<-rbind(matr, NA)
   j <- 1:ncol(matr)
   segments(intervals[rango], matr[rango,j], intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2), matr[rango,j], col=color, ...)  #### line with the y level
   segments(intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2), matr[rango,j], intervals[rango+1], matr[rango+1,j], col=color, lty="dotted", ...)  #### line joining the y levels
   abline(v=middle[rango], lty="dotted", col=colorlines)
   mtext(textAxis, side=2, adj=0.5, line=1.5, outer=FALSE, ...)
   axis(1, at=middle[1:length(rt)], labels=featureData(ecs)$exonID[rt], ...)
}


#########################
#FUNCTION TO WRITE THE GENE MODELS:
#########################
drawGene <- function(minx, maxx, tr, rango, exoncol=NULL, names, trName, ...)
{
   plot.new()
   plot.window(xlim=c(minx, maxx), ylim=c(0, 1))
   rect(tr[rango,1], 0, tr[rango,2], 1, col=exoncol)
   zr <- apply(rbind(tr[rango, 2], tr[rango+1, 1]), 2, median)
   segments(tr[rango,2], 0.5, zr, 0.65)
   segments(zr, 0.65, tr[rango+1,1], 0.5)
   if(names){
      mtext(trName, side=2, adj=0.5, padj=1, line=1, outer=FALSE, las=2, ...)
   }
}

