rna-seq-diff-exprn
=============

RNA-sequencing differential expression analysis pipeline.
---------------------------------------------------------

Unlike DNA, RNA cannot be sequenced directly because it is less stable.
RNA-sequencing (Nagalakshmi et al, Science, 2008, 
doi:10.1126/science.1158441) is a method of the unstable ribonucleic acid 
(RNA) into "coding-DNA" (cDNA) products that can be sequenced by a 
sequencing machine such as Illumina Hi-Seq 2000.

Performs: genome coverage (via bedtools and HTSeq), generates Circos code and plots, differential expression (via DESeq and NOISeq), structural variant detection (e.g. fusion genes, via SVDetect) and differential exon usage (via DEXSeq).

To run the example, go to the folder for `rna-seq-diff-exprn`. In my case, this is `/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/`. Then, perform this command:
```
scripts/pipeline.sh test-results test-data/conditions_chr9.tab test-data/hg19_ucsc_genes.gtf test-data/hg19_ucsc_genes_chr9_dexseq.gtf test-data/hg19_ucsc_genes.bed test-data/hg19_id_symbol.txt test-data/human.hg19.genome test-data/karyotype/karyotype.human.hg19.txt test-data/hg19_gene_density_1e5bins.txt test-data/hg19_gc_content_circos_chr9.txt 2
```

Dependencies (what you should already have installed, or need to install to use this RNA-Sequencing analysis software)
Note: This may seem like a lot, but if you are in biomedical research, it 
is likely that the servers at your institution already have most of these
installed.

1. Python 2.7 (for RSeQC and HTSeq, #2 and #3)
   http://www.python.org/getit/releases/2.7/
   1. Cython - required for pybedtools below
      http://www.cython.org/#download
   2. pybedtools - required for HTSeq
      http://packages.python.org/pybedtools/
   3. HTSeq - Another method of calculating gene expression counts,
      in addition to BEDTools Coverage.
      http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
      Installed via python2.7, i.e. instead of typing
       `python setup.py install`
     Say:
       `python2.7 setup.py install`

2. RSeQC, RNA-Sequencing Quality Control software 
   Version 2.3 or later, http://code.google.com/p/rseqc/downloads/list

3. BEDTools, version 2.16.2 or greater
   http://code.google.com/p/bedtools/downloads/list

4. Perl 5.8.x or newer (for Circos, #6)
   Linux and Mac users: http://www.perl.org/get.html
   Windows users: http://strawberryperl.com/

5. Circos, version 0.60 or later (for plotting genome coverage data)
   Note: Circos has a number of Perl package dependencies that take some
   time to make sure they are all properly installed on your system.
   `http://circos.ca/software/download/circos`
   Aliased such that `circos` will run the program
   On my machine, this is accomplished by adding this line to the file in
   `/Users/olgabotvinnik/.bashrc`, or my `~/.bashrc` file:
     `PATH=$PATH:/usr/bin/circos/bin ; export PATH`
   Which means that when you run commands, the computer will know to look 
   in `/usr/bin/circos/bin` for potential executable files. /usr/bin/circos
   Is where I personally installed Circos. On the server that I use, for
   example, it is installed in:
     `/share/apps/circos-0.60/bin`
   So then my `~/.bashrc` file on the server looks like:
     `PATH=$PATH:/share/apps/circos-0.60/bin ; export PATH`

6. R, 2.14.2 or later
   http://www.r-project.org/
   1. DESeq
   	  For differential expression analysis
      http://www-huber.embl.de/users/anders/DESeq/
   2. NOISeq
      Another type of differential expression analysis software
      http://bioinfo.cipf.es/noiseq/doku.php?id=downloads
   3. DEXSeq
      Differential exon usage analysis
      http://bioconductor.org/packages/release/bioc/html/DEXSeq.html


Other useful functions:
* `wig_to_circos.R`: Converts .wig files (genome browser-type files) to
  files compatible with the Cirocs graphing format
* `get_gene_density.R`: Using a knownCanonical format file, which looks like:

```
    #chrom  chromStart  chromEnd  clusterId transcript  protein
    chr1  11873 14409 1 uc010nxq.1  uc010nxq.1
    chr1  14361 19759 2 uc009viu.3  uc009viu.3
    chr1  14406 29370 3 uc009viw.2  uc009viw.2
    chr1  34610 36081 4 uc001aak.3  uc001aak.3
    chr1  69090 70008 5 uc001aal.1  uc001aal.1
    chr1  136697  140566  6 uc001aam.4  uc001aam.4
    chr1  321083  321115  7 uc001aaq.2  uc001aaq.2
    chr1  321145  321207  8 uc001aar.2  uc001aar.2
    chr1  322036  326938  9 uc009vjk.2  uc009vjk.2
  Find the number of genes per 1 megabase and reports the gene density, 
  for example:
    #chr  chrStart  chrEnd  density
    chr1 11873 1011872 0.00062000062000062
    chr1 1011873 2011872 0.000585000585000585
    chr1 2011873 3011872 0.00058000058000058
    chr1 3011873 4011872 0.000592000592000592
    chr1 4011873 5011872 0.000394000394000394
    chr1 5011873 6011872 0.0002000002000002
    chr1 6011873 7011872 0.00025000025000025
    chr1 7011873 8011872 0.000254000254000254
    chr1 8011873 9011872 0.000179000179000179
```