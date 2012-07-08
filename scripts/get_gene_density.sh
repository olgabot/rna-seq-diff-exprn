#!/bin/sh

# get_gene_density.sh
# Created by: Olga Botvinnik on 07-07-2012
# Purpose of this script: Given a list of genes and their locations,
#                         find the gene density at each 1 Mb chromosome 
#						  interval
# Input(s): knownCanonical gene file of this format:
# #chrom	chromStart	chromEnd	clusterId	transcript	protein
# chr1	11873	14409	1	uc010nxq.1	uc010nxq.1
# chr1	14361	19759	2	uc009viu.3	uc009viu.3
# chr1	14406	29370	3	uc009viw.2	uc009viw.2
# chr1	34610	36081	4	uc001aak.3	uc001aak.3
# chr1	69090	70008	5	uc001aal.1	uc001aal.1
# chr1	136697	140566	6	uc001aam.4	uc001aam.4
# chr1	321083	321115	7	uc001aaq.2	uc001aaq.2
# chr1	321145	321207	8	uc001aar.2	uc001aar.2
# chr1	322036	326938	9	uc009vjk.2	uc009vjk.2
#
# Output(s): Circos-formatted gene density file like this:
# chr1	11873	1011873	0.5
# ...

while read line ; do
	# remove everything with underscores and chrM
	THIS_CHR=

done <$1 >$1.density.txt