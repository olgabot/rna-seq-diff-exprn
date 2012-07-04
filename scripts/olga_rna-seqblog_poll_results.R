pdf("~/workspace/UCSC/pourmand/2012-Spring/writeup/figures/olga_rna-seqblog_poll_results.pdf")
pie(c(89, 11), labels=c("89% No", "11% Yes"), col=brewer.pal(3, "Set1")[1:2], #main="Do we yet have a firm handle on the most appropriate/accurate\n pipeline for analysis of RNA-Seq datasets?", 
sub="(n=72) source: rna-seqblog.com reader poll"
)
dev.off()

pdf("~/workspace/UCSC/pourmand/2012-Spring/writeup/figures/olga_rna-seqblog_poll_reasons.pdf", width=11, height=8.5)
par(las=2) # make label text perpendicular to axis
par(mar=c(5,25,4,2)) # increase y-axis margin.
reasons = c(38, 10, 10, 10, 8, 8, 7, 5, 4)
names(reasons) = c("Too many tools/methods - don't know where to begin",
	"Assembling short reads into full length transcripts",
	"Not sure how to interpret results",
	"Dealing with alternative splicing events",
	"Amount of data is overwhelming - difficult to handle",
	"Estimating transcript abundance",
	"Selecting appropriate normalization method",
	"Dealing with sequence bias in the data",
	"Lack of genome information for transcript mapping")
barplot(reasons, 
#	main="What is the biggest challenge when performing\nRNA-Seq data analysis?", 
	horiz=TRUE, 
	names.arg=paste(names(reasons), " (", reasons, "%)", sep=""),
	col=brewer.pal(3, "Dark2")[1], sub="(n=100) source: rna-seqblog.com reader poll")
dev.off()