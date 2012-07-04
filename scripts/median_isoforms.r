
args = commandArgs(TRUE)
combined.counts = args[1]

#counts = read.delim("ref_combined_counts.txt")
counts = read.delim(combined.counts)
genes.unique = unique(counts[,2])
counts.symbol = sapply(genes.unique, function(g) {
  ind = which(counts[,2]==g)
  c(median(counts[ind,3]),
    median(counts[ind,4]),
    var(counts[ind,3]),
    var(counts[ind,4]) ) } )
