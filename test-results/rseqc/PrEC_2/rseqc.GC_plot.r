pdf("/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/rseqc/PrEC_2/rseqc.GC_plot.pdf")
gc=rep(c(76.67,86.67,63.33,10.00,96.67,26.67,0.00,36.67,66.67,80.00,93.33,43.33,70.00,13.33,50.00,53.33,90.00,73.33,30.00,20.00,83.33,60.00,40.00,46.67,23.33,56.67,100.00,16.67,33.33),times=c(400,34,2543,4,8,524,1,1977,1912,169,4,3138,1294,15,3816,3840,12,740,951,133,79,3173,2421,3553,276,3683,4,51,1367))
hist(gc,probability=T,breaks=100,xlab="GC content (%)",ylab="Density of Reads",border="blue",main="")
dev.off()
