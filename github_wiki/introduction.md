I have worked on this package for the past 4 months and have put a lot of effort into making it as easy to use as possible, so please let me know if you run into any issues.

My vision of how easy it should be to use is outlined in the headline photo and title: 
> `.bam` files —-> ??? —> Results!

Structural variant detection is grey because it is not yet available.

The example included in the code looks at the differential expression of two prostate cancer cell lines, LNCaP and PrEC. It only includes reads from chromosome 9, an arbitrary choice to reduce the time of running the example. Currently, the example takes 20 minutes to run on my laptop (MacBook Pro, Late 2008 model, 2.4G Hz Intel Core Duo, 4 GB RAM), and this is really supposed to be run on a server that holds all your gigantic (10gb+) `.bam` files. If you do choose to run it on your own personal computer, I would recommend leaving your computer alone while it runs, or else your entire machine may crash. Mine has.


# The inputs

Basically, the inputs are (they will be more thoroughly described below):

1. Output directory
2. Metadata file with the `.bam` file locations, sample IDs, and experiment descriptions
3. Gene transfer format file (GTF file, a subset of the General Feature Format, GFF filetypes) of genes for your species
4. GTF file processed with DEXSeq, a differential exon usage program.
5. Browser-extensible data (BED) file of genes for your species. Both the BED and GTF files are needed for the different transcript-counting programs (`bedtools Coverage` and `htseq-count`)
6. Transcript ID - Symbol file (explained below)
7. "Genome" file. This is available in your BEDTools distribution for mouse and human (mm8, mm9, hg18 and hg19), and I'll show you an example of one so you can create one for other species, too.
8. "Karyotype" file for [Circos](http://circos.ca/) plotting. This is available in your Circos distribution for mouse and human (mm8, mm9, hg18 and hg19), and again, I'll show you an example so you can create your own.
9. Gene density file, for [Circos](http://circos.ca/) plotting. I have created one for humans (hg19) and will give instructions and a script for how to create one on your own.
10. GC content file, for [Circos](http://circos.ca/) plotting. I have one for humans, chromosome 9 (it takes a long time to plot if it is very dense, so I chose a subset of the genome), and will show you how to create your own
11. Number of "bunches" you'd like to create (more on this later)

# The outputs

And after submitting all those files and waiting ~20 min for the example to finish, you get:

1. Quality-control plots of the `.bam` files, via [RSeQC](http://code.google.com/p/rseqc/). Most of the commands listed in the [manual](http://code.google.com/p/rseqc/wiki/Manual) are performed. In alphabetical order, they are (with documentation copy/pasted from the RSeQC manual):
    * [`clipping_profile.py`](http://code.google.com/p/rseqc/wiki/Manual#clipping_profile.py): This program is used estimate clipping profile of RNA-seq reads from BAM or SAM file. Note that to use this function, CIGAR strings within SAM/BAM file should have `S` operation (This means your reads mapper should support clipped mapping).
    * [`geneBody_coverage.py`](http://code.google.com/p/rseqc/wiki/Manual#geneBody_coverage.py): Read coverage over gene body. This module is used to check if reads coverage is uniform and if there is any 5’/3’ bias. This module scales all transcripts to 100 nt and calculates the number of reads covering each nucleotide position. Finally, it generates a plot illustrating the coverage profile along the gene body. NOTE: this module requires lots of memory for large BAM files, because it load the entire BAM file into memory. 
        * This is a VERY useful function to look at how fully covered the transcripts are in your experiment. I highly recommend taking a look at the gene body coverage before interpreting differential expression.
    * [`infer_experiment.py`](http://code.google.com/p/rseqc/wiki/Manual#geneBody_coverage.py): This program is used to speculate how RNA-seq sequencing were configured, especially how reads were stranded for strand-specific RNA-seq data, through comparing reads' mapping information to the underneath gene model. 
    	1. Note: currently, RPKM_count.py does not take in the infer_experiment data but will in the future.
    * [`inner_distance.py`](http://code.google.com/p/rseqc/wiki/Manual#inner_distance.py): This module is used to calculate the inner distance (or insert size) between two paired RNA reads. The distance is the mRNA length between two paired fragments.
    * [`junction_annotation.py`](http://code.google.com/p/rseqc/wiki/Manual#junction_annotation.py): For a given alignment file (`-i`) in BAM or SAM format and a reference gene model (`-r`) in BED format, this program will compare detected splice junction to reference gene model.
    * [`junction_saturation.py`](http://code.google.com/p/rseqc/wiki/Manual#junction_saturation.py): It's very important to check if current sequencing depth is deep enough to perform alternative splicing analyses. For a well annotated organism, the number of expressed genes in particular tissue is almost fixed so the number of splice junctions is also fixed,　all splice junctions can be predetermined according to reference gene model. All (annotated) splice junctions should be rediscovered from a saturated RNA-seq data, otherwise, downstream alternative splicing analysis is problematic because low abundance splice junctions are missing. This module checks for saturation by resampling 5%, 10%, 15%, ..., 95% of total alignments from BAM or SAM file, and then detects splice junctions from each subset and compares them to reference gene model.
    * [`read_distribution.py`](http://code.google.com/p/rseqc/wiki/Manual#read_distribution.py): Provided a BAM/SAM file and reference gene model, this module will calculate how mapped reads were distributed over genome feature (like CDS exon, 5'UTR exon, 3' UTR exon, Intron, Intergenic regions). When genome features are overlapped (e.g. a region could be annotated as both exon and intron by two different transcripts) , they are prioritize as: CDS exons > UTR exons > Introns > Intergenic regions, for example,if a read was mapped to both CDS exon and intron, it will be assigned to CDS exons.
    * [`read_duplication.py`](http://code.google.com/p/rseqc/wiki/Manual#read_duplication.py): Two strategies were used to determine reads duplication rate:
        * Sequence based: reads with exactly the same sequence content are regarded as duplicated reads.
        * Mapping based: reads mapped to the same genomic location are regarded as duplicated reads. For splice reads, reads mapped to the same starting position and splice the same way are regarded as duplicated reads.
    * [`read_GC.py`](http://code.google.com/p/rseqc/wiki/Manual#read_GC.py)
        * Finds the GC content of your reads and outputs a histogram plot.
    * [`read_NVC.py`](http://code.google.com/p/rseqc/wiki/Manual#read_NVC.py): This module is used to check the nucleotide composition bias. Due to random priming, certain patterns are over represented at the beginning (5'end) of reads. This bias could be easily examined by NVC (Nucleotide versus cycle) plot. NVC plot is generated by overlaying all reads together, then calculating nucleotide composition for each position of read (or each sequencing cycle). In ideal condition (genome is random and RNA-seq reads is randomly sampled from genome), we expect A%=C%=G%=T%=25% at each position of reads.
        * I also highly recommend looking at this plot and possibly trimming the first 6-9bp of your reads.
    * [`read_quality.py`](http://code.google.com/p/rseqc/wiki/Manual#read_quality.py): According to SAM specification, if Q is the character to represent "base calling quality" in SAM file, then Phred Quality Score = `ord(Q) - 33`. Here `ord()` is python function that returns an integer representing the Unicode code point of the character when the argument is a unicode object, for example, `ord('a')` returns 97. Phred quality score is widely used to measure "reliability" of base-calling, for example, phred quality score of 20 means there is 1/100 chance that the base-calling is wrong, phred quality score of 30 means there is 1/1000 chance that the base-calling is wrong. In general: Phred quality score = `-10*log10P`, here `P` is probability that base-calling is wrong.
    * [`RPKM_count.py`](http://code.google.com/p/rseqc/wiki/Manual#RPKM_count.py): Given a BAM file and reference gene model, this program will calculate the raw count and RPKM values for transcript at exon, intron and mRNA level. For strand specific RNA-seq data, program will assign read to its parental gene according to strand rule, if you don't know the strand rule, run `infer_experiment.py`. Please note that chromosome ID, genome cooridinates should be concordant between BAM and BED files.
    * [`RPKM_saturation.py`](http://code.google.com/p/rseqc/wiki/Manual#RPKM_saturation.py): The precision of any sample statitics (RPKM) is affected by sample size (sequencing depth); “resampling” or “jackknifing” is a method to estimate the precision of sample statistics by using subsets of available data. This module will resample a series of subsets from total RNA reads and then calculate RPKM value using each subset. By doing this we are able to check if the current sequencing depth was saturated or not (or if the RPKM values were stable or not) in terms of genes' expression estimation. If sequencing depth was saturated, the estimated RPKM value will be stationary or reproducible. By default, this module will calculate 20 RPKM values (using 5%, 10%, ... , 95%,100% of total reads) for each transcripts. Although people can use KPSS test to determine if the estimated RPKM level is in stationary or not. Visual inspection is more accurate.
2. Genome-wide coverage of all samples, visualized by a [Circos](http://circos.ca) plot. The example only includes reads mapped to chromosome 9, so the example plot only shows coverage from chromosome 9.
3. Transcript counts, counted by two methods:
    * BEDTools, via `bedtools coverage`
    * HTSeq, via `htseq-count`
4. Differential Expression analysis, by `DESeq` (in the `R` programming language). This include gene lists and heat maps.
    * Other differential expression methods (NOISeq, DEXSeq) are in the works
5. (Not yet: Structural Variant detection, e.g. of fusion genes via SVDetect)