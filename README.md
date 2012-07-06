rna-seq-diff-exprn
====

RNA-sequencing differential expression analysis pipeline.

Unlike DNA, RNA cannot be sequenced directly because it is less stable.
RNA-sequencing (Nagalakshmi et al, Science, 2008, 
doi:10.1126/science.1158441) is a method of the unstable ribonucleic acid 
(RNA) into "coding-DNA" (cDNA) products that can be sequenced by a 
sequencing machine such as Illumina Hi-Seq 2000.

Performs: genome coverage (via bedtools and HTSeq), generates Circos code and plots, differential expression (via DESeq and NOISeq), structural variant detection (e.g. fusion genes, via SVDetect) and differential exon usage (via DEXSeq).

Dependencies (what you should already have installed, or need to install to use this RNA-Sequencing analysis software)
Note: This may seem like a lot, but if you are in biomedical research, it 
is likely that the servers at your institution already have most of these
installed.
1. Python 2.7 (for RSeQC and HTSeq, #2 and #3)
   http://www.python.org/getit/releases/2.7/
2. RSeQC, RNA-Sequencing Quality Control software 
   Version 2.3 or later, http://code.google.com/p/rseqc/downloads/list
3. HTSeq, installed via python2.7, i.e. instead of typing
     python setup.py install
   Say:
     python2.7 setup.py install
   http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
4. BEDTools, version 2.16.2 or greater
   http://code.google.com/p/bedtools/downloads/list
5. Perl 5.8.x or newer (for Circos, #6)
   Linux and Mac users: http://www.perl.org/get.html
   Windows users: http://strawberryperl.com/
6. Circos, version 0.60 or later (for plotting genome coverage data)
   Note: Circos has a number of Perl package dependencies that take some
   time to make sure they are all properly installed on your system.
   http://circos.ca/software/download/circos
   Aliased such that 'circos' will run the program
   On my machine, this is accomplished by adding this line to the file in
   /Users/olgabotvinnik/.bashrc, or my ~/.bashrc file:
     PATH=$PATH:/usr/bin/circos/bin ; export PATH
   Which means that when you run commands, the computer will know to look 
   in /usr/bin/circos/bin for potential executable files. /usr/bin/circos
   Is where I personally installed circos. On the server that I use, for
   example, it is installed in:
     /share/apps/circos-0.60/bin
   So then my ~/.bashrc file would look like:
     PATH=$PATH:/share/apps/circos-0.60/bin ; export PATH
7. R, 2.14.2 or later
   http://www.r-project.org/
   a. DESeq
   	  For differential expression analysis
      http://www-huber.embl.de/users/anders/DESeq/
   b. NOISeq
      Another type of differential expression analysis software
      http://bioinfo.cipf.es/noiseq/doku.php?id=downloads
   c. DEXSeq
      Differential exon usage analysis
      http://bioconductor.org/packages/release/bioc/html/DEXSeq.html