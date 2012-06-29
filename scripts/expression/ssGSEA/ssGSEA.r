#!/share/apps/bin/R

source("make.pinkogram.R")

############# Import arguments
args = commandArgs(TRUE)
results/expression.table = read.delim(args[1], header=TRUE)
nFeatures = 20
##############################

gene.names = make.unique(results/expression.table[,1])
results/expression.counts = results/expression.table[,-1]

mycol = make.pinkogram()
n.colors = length(mycol)
max.c <- max(max(results/expression.counts), -min(results/expression.counts))
heatm <- ceiling(n.colors * (results/expression.counts
                             - (- max.c))/(1.001*(max.c - (- max.c))))
heatm = heatm[1:nFeatures,]

nSample = dim(heatm)[2]
pdf("/home/obot/single-cell/results/expression/results/expression.pdf")
image(1:nSample, 1:nFeatures, heatm, zlim = c(0, n.colors),
      col=mycol, axes=FALSE, main ="", sub="", xlab="", ylab="")
axis(2, at=nFeatures, labels=gene.names[1:nFeatures], adj= 0.5,
     tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
axis(1, at=1:nSample, labels=col.names(heatm), adj= 0.5, tick=FALSE, las = 3,
     cex.axis=0.60, font.axis=1, line=-1)

# May need this:
## Tissue Legend
ifMultipleTissues = FALSE
if(ifMultipleTissues){
  #browser()
  par(mar = c(0, 0, 0, 0))
  plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
  legend(x=0, y= 1, #horiz = T, x.intersp = 0.5, y.intersp=.25,
         legend=tissue.names, bty="n", xjust=0, yjust= 1,
         fill = tissue.colors, cex = 0.8, #pt.cex=1.75,
         ncol=4)
}
