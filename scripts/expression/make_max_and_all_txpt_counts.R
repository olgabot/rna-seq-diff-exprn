#!/usr/bin/Rscript
options(error=utils::recover)


##
## 
# Takes in a dataset and ID symbols and outputs two datasets:
# 1. (MAX) For each gene symbol, take the max
# 2. (ALL) For each gene symbol, "unique-ify", 
#       ie many TP53s --> TP53, TP53.1, TP53.2, etc
# INPUT: dataframe of read counts per gene
#        dataframe with transcript ID and symbol
# OUTPUT: two datasets:
#		one that has the max for each gene symbol
#		one that has all possible transcripts for that gene symbol
##
## Created by  on .
##
## Collaborators: 
## Usage: ./.R <arg1> <arg2> <arg3>
## Example run: ./make_max_and_all_txpt_counts.R ../../test-results/expression/bedtools/bedtools_counts_table.tab ../../test-results/expression/bedtools/bedtools_counts_table_all.tab ../../test-results/expression/bedtools/bedtools_counts_table_max.tab


###### Get initial conditions ######
args = commandArgs(TRUE)
# testing:
# counts.ds = read.delim("~/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/bedtools_counts_table.tab", stringsAsFactors=FALSE)

counts.dsFile = args[1]
# idToSymbol = read.delim(args[2])

# The output filename of the "all" file, which takes all the 
# transcripts for a gene symbol and renames them to the gene symbol
# For example, if we start with these counts of the gene TTLL11:
# uc004blr.3	TTLL11	397	0	456	450	298	207	1744	838
# uc004blt.1	TTLL11	78	0	92	98	76	53	356	166
# uc004blu.1	TTLL11	78	0	92	98	76	53	356	166
# uc011lyl.2	TTLL11	397	0	456	450	298	207	1744	838
# uc011lym.1	TTLL11	71	0	88	87	61	37	322	147
# Then we end with these counts:
# TTLL11 (uc004blr.3)	397	0	456	450	298	207	1744	838
# TTLL11 (uc004blt.1)	78	0	92	98	76	53	356	166
# TTLL11 (uc004blu.1)	78	0	92	98	76	53	356	166
# TTLL11 (uc011lyl.2)	397	0	456	450	298	207	1744	838
# TTLL11 (uc011lym.1)	71	0	88	87	61	37	322	147
#  This way, when you look at the top differentially expressed
# transcripts, you know which transcript to look up.
countsAllFile = args[2]

# The output filename of the "max" file, which takes the maximum 
# counts for all the transcripts of a gene symbol for each sample.
# For example, if we start with these counts:
# uc004blr.3	TTLL11	397	0	456	450	298	207	1744	838
# uc004blt.1	TTLL11	78	0	92	98	76	53	356	166
# uc004blu.1	TTLL11	78	0	92	98	76	53	356	166
# uc011lyl.2	TTLL11	397	0	456	450	298	207	1744	838
# uc011lym.1	TTLL11	71	0	88	87	61	37	322	147
# Then we get these counts at the end:
# TTLL11	397	0	456	450	298	207	1744	838
# Which was the maximum for each sample, so we can't associate it
# with a particular transcript ID.
countsMaxFile = args[3]
####################################
print(ls.str())

counts.ds = read.delim(counts.dsFile, stringsAsFactors=FALSE)

geneSymbolsMax = unique(counts.ds[,2])

max.ds = t(sapply(geneSymbolsMax, function(symb) {
	symbInd = which(counts.ds[,2] == symb)

	# Cut out the first and second columns because those are
	# the gene names and we only want to take the max of the counts
	symbCounts = apply(counts.ds[symbInd,-1:-2], MARGIN=2, FUN=max)
	}))

geneSymbolsAll = paste(counts.ds[,2], " (", 
	counts.ds[,1], ")", sep="")
all.ds = counts.ds[,-1]
all.ds[,1] = geneSymbolsAll

header = paste(c("Name", colnames(counts.ds)[-1:-2]), collapse="\t")

write(header, file=countsAllFile, append=FALSE)
# --- row.names = FALSE for "all" because we made --- #
# --- the first column contain the gene symbols   --- #
write.table(all.ds, file=countsAllFile, append=TRUE, sep="\t",
	quote = FALSE, eol = "\n", na = "NA", dec = ".", 
	row.names = FALSE, col.names = FALSE)

write(header, file=countsMaxFile, append=FALSE)
# --- row.names = TRUE for "max" because the we don't keep --- #
# --- the gene symbol in the first column like in "all"    --- #
write.table(max.ds, file=countsMaxFile, append=TRUE, sep="\t",
	quote = FALSE, eol = "\n", na = "NA", dec = ".", 
	row.names = TRUE, col.names = FALSE)