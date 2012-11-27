# Input files: A thorough description

## Output directory

This is where all the results files will be. Folders that will be created are:

The basic directories are:
```
expression/
rseqc/
circos/
```

Quality control files go in:
```
rseqc/[all your sample ids]
```

Circos files go in:
```
circos/[all your sample ids]
```

The counts files will go in:
```
expression/bedtools/[all your sample ids]
expression/htseq/[all your sample ids]
```

Heatmaps will go:
```
expression/bedtools/figures/
expression/hsteq/figures/
```

## Metadata file

* What: Description of the input data files, where they are, sample IDs, a "group" such as a treatment type, and details of the RNA-Seq experiment
* Used by: Entire pipeline
* Example file: `test-data/conditions_chr9.tab`

Text of example file:
```
# These data are a subset of http://0-www.ncbi.nlm.nih.gov.elis.tmu.edu.tw/geo/query/acc.cgi?acc=GSE27619 from GEO, specifically:
# GSM721116, GSM721117, GSM721118, GSM721119, GSM721123, GSM721124
## LNCaP and PrEC are the names of prostate cancer cell lines used in this GEO dataset
## As you may have noticed, any lines starting with a hash (`#') are ignored
## This file is based on the columns, so if you add extraneous columns, they will not be read into the program (ie you can add comments)
####
# This is a tab-delimited file, with the following columns:
# bam_prefix     id     group     gender     read_type     strandedness
## bam_prefix = the filename, without the .bam
## id = (user-defined, could be anything)
## group = (user-defined, could be anything)
## gender = male --or-- female
## read_type = single_read --or-- paired_end
## strandedness = whether the cDNA library preparation was strand-specific or not. values = "strand_specific" --or-- "not_strand_specific"
test-data/GSM721117_mctp_20F0GAAXX_1_chr9_withheader     LNCaP_1     LNCaP     male     single_read     not_strand_specific # (This is a comment) Illumina
test-data/GSM721119_mctp_20F0GAAXX_2_chr9_withheader     LNCaP_2     LNCaP     male     single_read     not_strand_specific     # Illumina
test-data/GSM721118_mctp_20F0GAAXX_3_chr9_withheader     LNCaP_3     LNCaP     male     single_read     not_strand_specific
test-data/GSM721116_mctp_20F0GAAXX_4_chr9_withheader     LNCaP_4     LNCaP     male     single_read     not_strand_specific
test-data/GSM721123_mctp_30CYNAAXX_5_chr9_withheader     PrEC_1     PrEC     male     single_read     not_strand_specific
test-data/GSM721124_mctp_209ENAAXX_8_chr9_withheader     PrEC_2     PrEC     male     single_read     not_strand_specific
```

A tab-delimited file of the conditions of each sample, including:

1. `.bam` file location (without the final `.bam`, this makes the processing much easier)
    * Note: your `.bam` files MUST be indexed (have a accompanying `.bai` files) for RSeQC
2. sample ID that you can assign to whatever you want
3. The group or treatment type. In the example, it is the identification of two cell lines, but you may have "untreated" and "treated," for example.
4. Gender of the sample. This is important for Circos plotting, and whether we need to include or omit the Y chromosome.
5. Read type, whether the sequencing was single-read or paired-end
6. Strandedness, whether the cDNA library preparation was strand-specific or not.

## Gene Transfer Format (GTF) file for estimating gene counts

* What: GTF (Gene Transfer Format) files that you want to use to estimate gene counts.
* Used by: `HTSeq` and `DEXSeq` (differential exon usage)
* Example file: `test-data/hg19_ucsc_genes.gtf`

First few lines of example file:
```
chr1     hg19_knownGene     exon     11874     12227     0.000000     +     .     gene_id "uc001aaa.3"; transcript_id "uc001aaa.3";
chr1     hg19_knownGene     exon     12613     12721     0.000000     +     .     gene_id "uc001aaa.3"; transcript_id "uc001aaa.3";
chr1     hg19_knownGene     exon     13221     14409     0.000000     +     .     gene_id "uc001aaa.3"; transcript_id "uc001aaa.3";
chr1     hg19_knownGene     exon     11874     12227     0.000000     +     .     gene_id "uc010nxr.1"; transcript_id "uc010nxr.1";
chr1     hg19_knownGene     exon     12646     12697     0.000000     +     .     gene_id "uc010nxr.1"; transcript_id "uc010nxr.1";
chr1     hg19_knownGene     exon     13221     14409     0.000000     +     .     gene_id "uc010nxr.1"; transcript_id "uc010nxr.1";
chr1     hg19_knownGene     start_codon     12190     12192     0.000000     +     .     gene_id "uc010nxq.1"; transcript_id "uc010nxq.1";
chr1     hg19_knownGene     CDS     12190     12227     0.000000     +     0     gene_id "uc010nxq.1"; transcript_id "uc010nxq.1";
chr1     hg19_knownGene     exon     11874     12227     0.000000     +     .     gene_id "uc010nxq.1"; transcript_id "uc010nxq.1";
chr1     hg19_knownGene     CDS     12595     12721     0.000000     +     1     gene_id "uc010nxq.1"; transcript_id "uc010nxq.1"; 
```

GTF files are a subset of GFF (General Feature Format) files. To get one of these files, follow the following steps: 
(there is probably a similar method to use the ENSEMBL website, but I am not familiar with it so I am giving these instructions that I myself have followed many times)

1. Go to http://genome.ucsc.edu/cgi-bin/hgTables
2. Choose your clade and organism of interest
3. Choose these settings:
    * group: "Genes and Gene Prediction Tracks"
    * track: (whatever you want, but these instructions are built on using "UCSC Genes." You can use Ensembl or other transcript IDs but then you will need to choose different columns from the kgXref file for the TranscriptID-Symbol file, so I recommend sticking with UCSC IDs for now)
    * table: "knownGene"
    * region: "genome"
    * output format: "GTF - gene transfer format"
4. output file: (whatever you want, but I suggest something informative like `hg19_ucsc_genes.gtf`)
    * Make sure to include the file extension (`.gtf`) in the filename
5. Press "get output"
    * A file will be downloaded to your "Downloads" folder that you will need to move to somewhere more permanent, such as where you keep your other gene information files.


## GTF file processed with DEXSeq

* What: GTF (Gene Transfer Format) file specially formatted for use with DEXSeq, which measures differential exon usage.
* Used by: DEXSeq 
    * [not yet fully running, but you need to specify this so the order of files doesn't get altered from what the program expects]
* Example file: `test-data/hg19_ucsc_genes_chr9_dexseq.gtf`

First few lines of example file:
```
chr9     hg19_ucsc_genes_chr9.gtf     aggregate_gene     11987     14525     .     +     .     gene_id "uc011llp.1"
chr9     hg19_ucsc_genes_chr9.gtf     exonic_part     11987     12340     .     +     .     transcripts "uc011llp.1"; exonic_part_number "001"; gene_id "uc011llp.1"
chr9     hg19_ucsc_genes_chr9.gtf     exonic_part     12726     12834     .     +     .     transcripts "uc011llp.1"; exonic_part_number "002"; gene_id "uc011llp.1"
chr9     hg19_ucsc_genes_chr9.gtf     exonic_part     13334     14525     .     +     .     transcripts "uc011llp.1"; exonic_part_number "003"; gene_id "uc011llp.1"
chr9     hg19_ucsc_genes_chr9.gtf     aggregate_gene     14511     29739     .     -     .     gene_id "uc010mgp.1+uc010mgm.1+uc022bcs.1+uc011llq.1+uc003zfu.1"
chr9     hg19_ucsc_genes_chr9.gtf     exonic_part     14511     14940     .     -     .     transcripts "uc010mgm.1"; exonic_part_number "001"; gene_id "uc010mgp.1+uc010mgm.1+uc022bcs.1+uc011llq.1+uc003zfu.1"
chr9     hg19_ucsc_genes_chr9.gtf     exonic_part     15081     15149     .     -     .     transcripts "uc010mgm.1+uc022bcs.1"; exonic_part_number "002"; gene_id "uc010mgp.1+uc010mgm.1+uc022bcs.1+uc011llq.1+uc003zfu.1"
chr9     hg19_ucsc_genes_chr9.gtf     exonic_part     15909     16061     .     -     .     transcripts "uc010mgm.1+uc022bcs.1"; exonic_part_number "003"; gene_id "uc010mgp.1+uc010mgm.1+uc022bcs.1+uc011llq.1+uc003zfu.1"
chr9     hg19_ucsc_genes_chr9.gtf     exonic_part     16188     16421     .     -     .     transcripts "uc011llq.1"; exonic_part_number "004"; gene_id "uc010mgp.1+uc010mgm.1+uc022bcs.1+uc011llq.1+uc003zfu.1"
chr9     hg19_ucsc_genes_chr9.gtf     exonic_part     16718     16876     .     -     .     transcripts "uc010mgm.1+uc022bcs.1+uc011llq.1"; exonic_part_number "005"; gene_id "uc010mgp.1+uc010mgm.1+uc022bcs.1+uc011llq.1+uc003zfu.1"
```

To create this file, You need a GTF file (which can be obtained as described in the GTF section), and then use the script included in `rna-seq-diff-exprn/scripts/external/dexseq_prepare_annotation.py`:
```
$ python2.7 scripts/external/dexseq_prepare_annotation.py test-data/hg19_ucsc_genes.gtf test-data/hg19_ucsc_genes_dexseq.gtf
```
The dollar sign `$` indicates a bash shell and shows that we are using a command-line interface, as opposed to a command embedded in source code such as this document.


## Browser extensible data (BED) files, again for estimating gene counts

* What: BED (Browser Extensible Data) files that you want to use to estimate gene counts.
* Used by: `bedtools coverage`
* Example file: `test-data/hg19_ucsc_genes.bed`

First few lines of example file:
```
track name="tb_knownGene" description="table browser query on knownGene" visibility=3 url=
chr1     11873     14409     uc001aaa.3     0     +     11873     11873     0     3     354,109,1189,     0,739,1347,
chr1     11873     14409     uc010nxr.1     0     +     11873     11873     0     3     354,52,1189,     0,772,1347,
chr1     11873     14409     uc010nxq.1     0     +     12189     13639     0     3     354,127,1007,     0,721,1529,
chr1     14361     16765     uc009vis.3     0     -     14361     14361     0     4     468,69,147,159,     0,608,1434,2245,
chr1     16857     17751     uc009vjc.1     0     -     16857     16857     0     2     198,519,     0,375,
chr1     15795     18061     uc009vjd.2     0     -     15795     15795     0     5     152,159,198,136,456,     0,811,1062,1437,1810,
chr1     14361     19759     uc009vit.3     0     -     14361     14361     0     9     468,69,152,159,198,510,147,99,847,     0,608,1434,2245,2496,2871,3553,3906,4551,
chr1     14361     19759     uc009viu.3     0     -     14361     14361     0     10     468,69,152,159,198,510,147,102,54,847,     0,608,1434,2245,2496,2871,3553,3906,4139,4551,
chr1     14361     19759     uc001aae.4     0     -     14361     14361     0     10     468,69,152,159,198,136,137,147,99,847,     0,608,1434,2245,2496,2871,3244,3553,3906,4551,
```

To get one of these files, do the following steps:
(there is probably a similar method to use the ENSEMBL website, but I am not familiar with it so I am giving these instructions that I myself have followed many times and can help you out with if you are stuck)

1. Go to http://genome.ucsc.edu/cgi-bin/hgTables
2. Choose your clade and organism of interest
3. Choose these settings:
    * group: "Genes and Gene Prediction Tracks"
    * track: (whatever you want, but these instructions are built on using "UCSC Genes." You can use Ensembl or other transcript IDs but then you will need to choose different columns from the kgXref file for the `TranscriptID-Symbol` file, so I recommend sticking with UCSC IDs for now)
    * table: "knownGene"
    * region: "genome"
    * output format: "BED - browser extensible data"
4. output file: (whatever you want, but I suggest something informative like `hg19_ucsc_genes.bed`)
    * Make sure to include the file extension (`.bed`) in the filename
5. Press "get output"
    * A file will be downloaded to your "Downloads" folder that you will need to move to somewhere more permanent, such as where you keep your other gene information files.

If you have an existing BED file, make sure the first line has `track name=….` or else RSeQC gets mad.


## TranscriptID-Symbol file

* What: Tab-delimited file with the transcript ID in column 1 and the gene symbols you'd like to use in column 2, without a header line.
* Used by: the pipeline to create a table of gene counts
* Example file: `test-data/hg19_id_symbol.txt`

First few lines of example file:
```
uc001aaa.3     DDX11L1
uc010nxr.1     DDX11L1
uc010nxq.1     DDX11L9
uc001aal.1     OR4F5
uc001aaq.2     DQ597235
uc001aar.2     DQ599768
uc001aau.3     LOC100132287
uc021oeh.1     LOC100133331
uc009vjk.2     LOC100133331
uc021oei.1     LOC388312
```

You can get one of these files by:

1. Go to http://genome.ucsc.edu/cgi-bin/hgTables
2. Choose your clade and organism of interest
3. Choose these settings:
    * group: "Genes and Gene Prediction Tracks"
    * track: (whatever you want, but these instructions are built on using "UCSC Genes." You can use Ensembl or other transcript IDs but then you will need to choose different columns from the kgXref file for the `TrascriptID-Symbol` file, so I recommend sticking with UCSC IDs for now)
    * table: "kgXref"
    * region: "genome"
    * output format: "GTF - gene transfer format"
4. output file: (whatever you want, but I suggest something informative like `hg19_kgXref.txt`)
    * Make sure to include the file extension (`.txt`, for example) in the filename
5. Press "get output"
6. Now you need to take an extra step to get just the UCSC IDs, e.g. `uc002gig.1` (column 1 in the kgXref file) and the gene symbols, e.g. `TP53` (column 5 in the kgXref file), a known tumor suppressor gene.

**But wait! You're not done yet!**

You need to remove the first line, the header of the file that explains what is in which column.

You could do this in Microsoft Excel, but the human file (for example) has 80,923 lines in it and will crash Excel. For organisms with fewer documented genes, using Excel to push columns around may be just fine.

The Linux/UNIX (lovingly called "\*NIX") commands to take columns is called "cut" (there is also "paste" to put together columns from different files but that's out of the scope of what we're doing here) We want columns 1 and 5 (the UCSC ID and the gene symbol - take a peek at the file by typing `head hg19_kgXref.txt` on the command line in the directory - this will show the first 10 lines of the file), so we'll say `cut -f 1,5` where the `-f` indicates the "fields" we want to `cut`. Then we use `sed 1d` to skip the first line (skipping more than one line has a slighly different command, check out [my favorite Sed tutorial](http://www.grymoire.com/Unix/Sed.html) if you're interested in learning more). And the `<` indicates the input file, the `|` indicates that the output of the previous command is treated as input to the next command, and the `>` indicates the output file. Note that we created a *new* file and did not overwrite the old one. In general, it is best practices to create a new file rather than overwrite the old one. Also, if you try to make your input and output files the same, the commands may get confused and you could lose your original data. :(
```
$ cut -f 1,5 < hg19_kgXref.txt | sed 1d > hg19_id_symbol.txt
```
Or, if you want to create a chromosome-specific file like I did, use your `.bed` file to search through your kgXref file:
```
$ cut -f 4 hg19_ucsc_genes_chr9.bed | grep --fixed-strings - hg19_kgXref.txt >hg19_kgXref_chr9.txt
```
Then do the same as above, but with your chr9 file:
```
$ cut -f 1,5 < hg19_kgXref_chr9.txt | sed 1d > hg19_id_symbol_chr9.txt
```


## Genome file

* What: This "Genome" file really just says how long each chromosome is.
* Used by: `genomeCoverageBed`
* Example file: `test-data/human.hg19.genome`

First few lines of example file:
```
chr1     249250621
chr2     243199373
chr3     198022430
chr4     191154276
chr5     180915260
chr6     171115067
chr7     159138663
chrX     155270560
chr8     146364022
chr9     141213431
```

Besides the example files, you can also use ones shipped with BEDTools. On my machine, these files are located in `~/packages/BEDTools/genomes`:
```
$ ls ~/packages/BEDTools-Version-2.16.2/genomes/
  human.hg18.genome human.hg19.genome mouse.mm8.genome  mouse.mm9.genome
```

As to how to create these files for non-mouse or human organisms, my suggestion (while rather unwieldy) is to:

1. Go to http://genome.ucsc.edu/cgi-bin/hgTables
2. Choose your clade and organism of interest
3. Choose these settings:
    * group: "All Tables"
    * table: "chromInfo"
    * output format: "all fields from selected table"
    * output file: (anything you want, but preferably something informative like platypus.ornAna1.genome)
4. Press "get output"
5. Remove the first line and the third column of the file, which you could do in Microsoft Excel (since this file will be comparatively small), or by using this shell command:
`$ cut -f 1,2 < platypus.ornAna1.genome | sed 1d > platypus.ornAna1.genome.fixed`


## Karyotype file

* What: Another "this is how long all the chromosomes are" file, but formatted so Circos can use it
* Used by: Circos
* Example file: `test-data/karyotype/karyotype.human.hg19.txt`

This file looks like:
```
chr - chr1 chr1 0 249250621 chr1
chr - chr2 chr2 0 243199373 chr2
chr - chr3 chr3 0 198022430 chr3
chr - chr4 chr4 0 191154276 chr4
chr - chr5 chr5 0 180915260 chr5
chr - chr6 chr6 0 171115067 chr6
chr - chr7 chr7 0 159138663 chr7
chr - chr8 chr8 0 146364022 chr8
chr - chr9 chr9 0 141213431 chr9
```

Karyotype file used by Circos, which specifies the chromosome lengths. The third column, the chromosome name, MUST use `chr1`-type notation, and not the typical `hs1` notation for Homo sapiens chromosome 1. This is because `bedtools` and friends use `chr1` notation, but I didn't want to require the organism name and then lookup the conversion. Presumably, you would have samples from all the same organism since you are comparing gene expression and coverage across different treatments, so I felt this was a safe assumption. I also didn't want to lock you into ONLY using human data, because there are plenty of interesting organisms out there.


## GC Content file

* What: Percentage of Guanine and Cytosine (GC) bases per some distance (I recommend at least 1000 bases, or else it will take a VERY long time to plot). This is an important quantity because the GC base pairing has three hydrogen bonds instead of two like AT, is known to be a stronger bond, and genes are also known to be GC-rich.
* Used by: Circos
* Example file: `test-data/hg19_gc_content_circos_chr9.txt`

First few lines of example file:
```
chr9     10000     10999     58
chr9     11000     11999     58
chr9     12000     12999     57.9
chr9     13000     13999     57.9
chr9     14000     14999     58
chr9     15000     15999     60.7
chr9     16000     16999     56.9
chr9     17000     17999     61.3
chr9     18000     18999     60.7
chr9     19000     19999     57.7
```

This GC content file can be created by converting a `.wig` (wiggle) format file that's used for a genome browser into a circus format using:
    ../scripts/wig_to_circos.R hg19_gc1000Base.txt hg19_gc_content.circos
    

## Number of bunches

* What: Merges samples of the same group, e.g. *untreated* into "bunches," or mergings of several samples. 
    * For example, in the example data provided, there are four samples of the LNCaP prostate cancer cell line, and two of the PrEC prostate cancer cell line.
    * In the provided example, two bunches are specified. This means that in addition to calculating the gene counts for `LNCaP_1`, `LNCaP_2`, `LNCaP_3`, `LNCaP_4`, `PrEC_1`, `PrEC_2`, the pipeline will also create 
    	* `LNCaP_bunch1of2`, containing `LNCaP_1` and `LNCaP_2`
    	* `LNCaP_bunch2of2`, containing `LNCaP_3` and `LNCaP_4`
    	* `PrEC_bunch1of2`, containing `PrEC_1`
    	* `PrEC_bunch2of2`, containing `PrEC_2`		
    * **Why this is useful:** Originally, this pipeline was written for single-cell RNA-Seq analysis and we found that while it is very interesting to look at the genes expressed in a single cell, we saw significantly more signal by merging treatment groups into bunches, and performing differential expression analysis on these bunches.
* Used by: Entire pipeline
* Example: 2

To change the ordering of samples, change the conditions file, as that is the basis of the bunching order.

If you don't want **any** bunches, omit this variable from the command line.