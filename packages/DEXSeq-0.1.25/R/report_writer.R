DEXSeqHTML <- function(ecs, geneIDs=NULL, path="DEXSeqReport", file="testForDEU.html", fitExpToVar="condition", FDR=0.1, color=NULL, color.samples=NULL, mart="", filter="", attributes="", extraCols=NULL)
{
   stopifnot(is(ecs, "ExonCountSet"))
   if(any(is.na(sizeFactors(ecs)))){
      stop("Please estimate sizeFactors first\n")}
   if(!fitExpToVar %in% ecs@designColumns){
      stop("fitExpToVar parameter is not in the design columns, double check ecs@designColumns")}
   if(sum(is.na(featureData(ecs)$dispersion))==nrow(counts(ecs))){
      stop("No dispersion parameters found, first call function fitDispersionFunction...\n")}

   ######## GET THE RESULT TABLE READY ##########
   results<-DEUresultTable(ecs)
   results[,c("dispersion", "pvalue", "padjust")] <- round(results[,c("dispersion", "pvalue", "padjust")], 3)
   if(!is.null(results$log2change)){
      results$log2change <- round(results$log2change, 3)
   }
   results[which(is.na(results$pvalue)),][,c("pvalue","padjust")]=1
   rownames(results) <- NULL
   numcond <- length(unique(design(ecs, drop=FALSE)[[fitExpToVar]]))
#	sortabletag <- hmakeTag(tag="script", src=(system.file(package="DEXSeq"), "/sorttable.js", sep=''))
   if(is.null(color)){
      color<-rgb(colorRamp(c("#D7191C", "#FFFFBF", "#2B83BA"))(seq(0, 1, length.out=numcond)), maxColorValue=255, alpha=175)
   }
   names(color) <- sort(levels(design(ecs, drop=FALSE)[[fitExpToVar]]))
   ### MAKING THE MAIN DIRECTORY AND files DIRECTORY
   dir.create(file.path(path, "files"), recursive=TRUE)
   ### SPECIFY THE GENES TO PUT IN THE REPORT ###

   if(is.null(geneIDs)){
      gns <- as.character(unique(results$geneID[which(results$padjust < FDR)]))
   }else{
      gns <- geneIDs
   }

   if(!all(gns %in% levels(geneIDs(ecs)))){
      stop("The geneIDs provided are not in the ecs object")}
   if(length(gns)==0){ 
      stop("There are no significant results in the test... nothing to report")}

   p<-openPage(file.path(path, file))
#   hwrite(sortabletag, p)
   hwrite('DEXSeq differential exon usage test', p, heading=1)
   hwrite('Experimental design', p, heading=2)
   cond<-as.matrix(cbind(rownames(design(ecs, drop=FALSE)), design(ecs, drop=FALSE)))
   rownames(cond) <- NULL
   colnames(cond)<-c("sample", ecs@designColumns)
#####
   condcolor <- matrix(rep("white", nrow(cond)*ncol(cond)), nrow(cond))
   condcolor[,which(colnames(cond) %in% fitExpToVar)] <- color[as.character(design(ecs, drop=FALSE)[[fitExpToVar]])]
#####
   if(!is.null(color.samples)){
      condcolor[,1] <- color.samples}
   hwrite(cond, bgcolor=condcolor, p)
   hwrite(paste("\n\nformulaDispersion = ", ecs@formulas[["formulaDispersion"]], sep=""), p, heading=3)
   hwrite(paste("\nformula0 = ", ecs@formulas[["formula0"]], sep=""), p, heading=3)
   hwrite(paste("\nformula1 = ", ecs@formulas[["formula1"]], sep=""), p, heading=3)
   hwrite('testForDEU result table', p, heading=2)
   ptowrite <- file.path(path, "files/")
   ######### prepare colors for table results
	
   m2col <- colorRampPalette(c("#FF706B", "#FEE08B", "white"), space="rgb")(5)
   matcol <- matrix(rep(m2col[5], ncol(results)*nrow(results)), nrow(results))
   ######### COLOR LEGEND = I BET THERE IS A SMARTER WAY OF DOING THIS:
   j <- 4
   for(i in c(0.25, 0.1, 0.05, 0.01)){
      matcol[which(results$padjust <= i),] <- m2col[j]
      j <- j-1
   }
   legend <- hwrite(c("<= 0.01", "<= 0.05", "<= 0.1", "<= 0.25", "> 0.25"), bgcolor=m2col)

   for( gene in gns){
      back <- hwrite("back", link=file.path("..", file))
      nameforlinks <- sapply(strsplit(gene, "\\+"), "[[", 1)
      otherlinks <- hwrite(c("counts", "expression", "splicing", "transcripts", "results"), link=c(paste(nameforlinks, "counts.html", sep=""), paste(nameforlinks, "expression.html", sep=""), paste(nameforlinks, "splicing.html", sep=""), paste(nameforlinks, "transcripts.html", sep=""), paste(nameforlinks, "results.html", sep="")), table=FALSE)
      loc <- as.character(results$geneID) %in% as.character(gene)
      ### this makes the page where to explore the pvalues ###
      subres <- results[loc,]
      submatcol <- matcol[loc,]
      rownames(subres) <- NULL
      genpage <- openPage(paste(ptowrite, nameforlinks, "results.html", sep=""))
      hwrite(c(back, otherlinks), table=TRUE, border=0, genpage)
      hwrite(legend, page=genpage)
      hwrite(as.matrix(subres), bgcolor=submatcol, table.class="sortable", style='margin:16px; border:0px solid black; border-width:0px; width:200px', table=TRUE, page=genpage)
      close(genpage, splash=TRUE)

		
      ### MAKE THE PLOT HTML PAGES FOR expression, counts and splicing
      makePlotPage(ecs=ecs, ptowrite=ptowrite, gene=gene, whichtag="expression", links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1200, height=700, h=7)
      makePlotPage(ecs=ecs, ptowrite=ptowrite, gene=gene, whichtag="counts", links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1200, height=700, h=7)
      makePlotPage(ecs=ecs, ptowrite=ptowrite, gene=gene, whichtag="splicing", links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1200, height=700, h=7)

      if(!is.null(featureData(ecs)$transcripts)){
         transcripts <- sapply(featureData(ecs)$transcripts[featureData(ecs)$geneID %in% gene], function(x){strsplit(x, ";")})
         trans <- Reduce(union, transcripts)
         h <- ifelse(length(trans) > 10, 7+(length(trans)*0.3), 7)
      if(sum(loc) > 30){   ############# if there are more than 30 exons, increase the size of the plotting region
         h <- h + (sum(loc)*0.1)
      }
      r <- try(
         makePlotPage(ecs=ecs, ptowrite=ptowrite, gene=gene, whichtag=c("expression", "transcripts"), links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1200, height=h*100, h=h)
           , silent=TRUE)
      }
   }
	
	
   results <- results[as.character(results$geneID) %in% gns,]
   genetable <- cbind( 
      geneID=unique(as.character(results$geneID)),
      do.call(rbind, 
         lapply(unique(as.character(results$geneID)), function(gene){
            vec <- as.character(geneIDs(ecs)) %in% gene 
            data.frame(chr=unique(featureData(ecs)$chr[vec]), start=min(featureData(ecs)$start[vec]), end=max(featureData(ecs)$end[vec]))
         })
      ),
      total_exons = rle(as.character(results$geneID))$lengths,
      exon_changes = sapply(unique(as.character(results$geneID)), function(gene){vec <- as.character(results$geneID) %in% gene; sum(results$padjust[vec] < FDR)})
   )

   if(class(mart) == "Mart"){
      if(attributes(mart)$dataset != ""){
      forvalues <- strsplit(as.character(genetable$geneID), "\\+")
      names(forvalues) <- genetable$geneID
      if(length(filter) > 1){
         warning("length(filter) > 2, only first element will be taken")
         filter <- filter[1]
      }
      extra <- getBM(attributes=c(filter, attributes), filters=filter, values=forvalues, mart=mart)
      fromart <- lapply(genetable$geneID, function(x){
         sep <- do.call(c, strsplit(as.character(x), "\\+"))
         extra[which(extra[,filter] %in% sep),]
      })

      extra <- sapply(attributes, 
         function(r){
           unlist(
              lapply(fromart, 
                 function(x){
                    paste(x[,r], collapse=" ")
                  }
              )
           )
         }
      )
      genetable <- cbind(geneID=genetable$geneID, extra, genetable[,2:length(genetable)])
      }else{
         warning("No dataset in biomart specified")
      }
   }else if( mart != ""){
      warning("Please provide a Mart class object for parameter mart")
   }  

   if( !is.null( extraCols )){
      genetable <- cbind( extraCols, genetable )
   }
   
	
   genetable$geneID <- sapply(as.character(genetable$geneID), function(m){w <- strsplit(m, "\\+");ns <- sapply(w, "[[", 1);hwrite(paste(unlist(w), collapse=" "), link=paste("files/", ns, "expression.html", sep=""))})
   rownames(genetable) <- NULL
   hwrite(genetable, page=p, table=TRUE, table.class="table-layout:fixed", style='margin:16px; border:0px solid black; border-width:1px; width:20%')
   close(p, splash=TRUE)
}

makePlotPage <- function(ecs, ptowrite, gene, whichtag, links, color, color.samples, FDR, fitExpToVar, width, height, h)
{
   allopts <- c("expression", "splicing", "counts", "transcripts")
   opts <- allopts %in% whichtag
   onlytag <- allopts[max(which(opts))]
   pagename <- sapply(strsplit(as.character(gene), "\\+"), "[[", 1)
   genpage <- openPage(paste(ptowrite, pagename, onlytag, ".html", sep=""))
   hwrite(links, table=TRUE, border=0, genpage)
   svg(paste(ptowrite, pagename, onlytag, ".svg", sep=""), height=h, width=12, pointsize=14)
   plotDEXSeq(ecs, geneID=gene, FDR=FDR, lwd=2, expression=opts[1], splicing=opts[2], norCounts=opts[3], displayTranscripts=opts[4], fitExpToVar=fitExpToVar, legend=TRUE, color=color, color.samples=color.samples, cex.axis=1.5)
   dev.off()
   hwrite(hmakeTag("iframe", src=paste(pagename, onlytag, ".svg", sep=""), width=width, height=height, border=0), page=genpage)
   close(genpage, splash=TRUE)
}
