# © 2012 Broad Institute, Inc.  All rights reserved.

#
# Projects gene expression data or ATARiS consistency scores (in either
# case, encoded in GCT input files) to gene set enrichment scores.
# (note: for clarity comments will mostly refer to expression data)
#
# Authors: Pablo Tamayo, Chet Birger
#

ssGSEA.project.dataset <- function(

    javaexec,

    jardir,

    # name of gct input file containing
    # gene expression data
    input.ds,

    # name of gct output file name containing
    # single-sample gene-set enrichment scores
    output.ds,

    # name of gene sets database from GSEA/MSigDB website
    gene.sets.database,

    # list of full pathnames of gmt input files
    # containing gene set definitions
    gene.sets.dbfile.list=NA,

    # column containing gene symbols within gct input file
    #  "Name", "Description"
    gene.symbol.column= "Name",

    # "ALL" or list with names of gene sets
    gene.set.selection  = "ALL",

    # normalization method applied to input feature data:
    # "rank", "log" or "log.rank"
    sample.norm.type    = "rank",

    # exponential weight applied to ranking in calculation of
    # enrichment score
    weight              = 0.75,

    # "combine.off" do not combine *_UP and *_DN versions in
    #   a single score. "combine.replace" combine *_UP and
    #   *_DN versions in a single score that replaces the individual
    # *_UP and *_DN versions. "combine.add" combine *_UP and
    # *_DN versions in a single score and add it but keeping
    # the individual *_UP and *_DN versions.
    combine.mode        = "combine.add",

    # min overlap required between genes in gene set and genes in input (feature
    # dataset) file in order to include that gene set in data set projection
    min.overlap         = 1)

{ #----------------------------------------------------------------------------------------

    # validate input parameters
    if (gene.symbol.column != "Name" && gene.symbol.column != "Description") {
        stop("invalid value for gene.symbol.column argument: ", gene.symbol.column)
    }
    if (sample.norm.type != "rank" && sample.norm.type != "log" && sample.norm.type != "log.rank") {
        stop("invalid value for sample.norm.type.argument: ", sample.norm.type)
    }
    if (combine.mode != "combine.off" && combine.mode != "combine.replace" && combine.mode != "combine.add") {
        stop("invalid value for combine.mode argument: ", combine.mode)
    }

    # Read input dataset (GCT format)
    dataset <- read.dataset(input.ds)
    m <- dataset$data

    # in Ataris or hairpin gct files the gene symbols are in the descs column
    if (gene.symbol.column == "Name") {
        gene.names <- row.names(m)
    }
    else if (gene.symbol.column == "Description") {
        gene.names <- dataset$row.descriptions
        row.names(m) <- gene.names
    }

    gene.descs <- dataset$row.descriptions
    sample.names <- colnames(m)

    Ns <- length(m[1,])
    Ng <- length(m[,1])

    # Sample normalization

    if (sample.norm.type == "rank") {
        for (j in 1:Ns) {  # column rank normalization
            m[,j] <- rank(m[,j], ties.method = "average")
        }
        m <- 10000*m/Ng
    } else if (sample.norm.type == "log.rank") {
        for (j in 1:Ns) {  # column rank normalization
            m[,j] <- rank(m[,j], ties.method = "average")
        }
        m <- log(10000*m/Ng + exp(1))
    } else if (sample.norm.type == "log") {
        m[m < 1] <- 1
        m <- log(m + exp(1))
    }

    # Read gene set databases

    if (!is.na(gene.sets.dbfile.list)) {

        # identify largest gene set size (max.G) and total number of
        # gene sets across all databases (max.N)
        max.G <- 0
        max.N <- 0
        for (gsdb in gene.sets.dbfile.list) {
            GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
            max.G <- max(max.G, max(GSDB$size.G))
            max.N <- max.N +  GSDB$N.gs
        }

        # create matrix (gs) containing all gene set definitions
        N.gs <- 0
        gs <- matrix("null", nrow=max.N, ncol=max.G)
        gs.names <- vector(length=max.N, mode="character")
        gs.descs <- vector(length=max.N, mode="character")
        size.G <- vector(length=max.N, mode="numeric")
        start <- 1
        for (gsdb in gene.sets.dbfile.list) {
            GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
            N.gs <- GSDB$N.gs
            gs.names[start:(start + N.gs - 1)] <- GSDB$gs.names
            gs.descs[start:(start + N.gs - 1)] <- GSDB$gs.desc
            size.G[start:(start + N.gs - 1)] <- GSDB$size.G
            gs[start:(start + N.gs - 1), 1:max(GSDB$size.G)] <- GSDB$gs[1:N.gs, 1:max(GSDB$size.G)]
            start <- start + N.gs
        }
        N.gs <- max.N
    }
    else {

        tmpGSDBFile <- tempfile()

        Get.GeneSets.db(javaexec, jardir, gene.sets.database, tmpGSDBFile)
        GSDB <- Read.GeneSets.db(tmpGSDBFile, thres.min=2, thres.max=2000, gene.names = NULL)
        unlink(tmpGSDBFile)

        max.G <- max(GSDB$size.G)
        N.gs <- GSDB$N.gs
        gs <- matrix("null", nrow=N.gs, ncol=max.G)
        size.G <- GSDB$size.G
        gs.names <- GSDB$gs.names
        gs.descs <- GSDB$gs.desc
        gs <- GSDB$gs[1:N.gs, 1:max.G]
    }

    # Select desired gene sets

    if (gene.set.selection[1] != "ALL") {
        locs <- match(gene.set.selection, gs.names)
        N.gs <- sum(!is.na(locs))
        if (N.gs == 0) {
            stop("No matches with gene.set.selection")
        }
        locs <- locs[!is.na(locs)]
        if(N.gs > 1) {
            gs <- gs[locs,]
        } else {
          # Force vector to matrix if only one gene set specified
            gs <- t(as.matrix(gs[locs,]))
        }
        gs.names <- gs.names[locs]
        gs.descs <- gs.descs[locs]
        size.G <- size.G[locs]
    }

    # Loop over gene sets

    # score.matrix records the enrichment score for each pairing
    # of gene set and sample
    score.matrix <- matrix(0, nrow=N.gs, ncol=Ns)
    for (gs.i in 1:N.gs) {
        gene.set <- gs[gs.i, 1:size.G[gs.i]]
        gene.overlap <- intersect(gene.set, gene.names)
        print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", length(gene.overlap)))
        if (length(gene.overlap) < min.overlap) {
            # if overlap between gene set and genes in input data set
            # are below a minimum overlap, no enrichment scores are
            # calculated for that gene set.
            score.matrix[gs.i, ] <- rep(NA, Ns)
            next
        }
        else {
            gene.set.locs <- match(gene.overlap, gene.set)
            gene.names.locs <- match(gene.overlap, gene.names)
            msig <- m[gene.names.locs,]
            msig.names <- gene.names[gene.names.locs]
            gs.score <- Project.to.GeneSet(data.array = m, weight = weight,
                                           gene.set = gene.overlap)
            score.matrix[gs.i,] <- gs.score$ES.vector

        }
    }

    # overlap pruning

    # eliminate gene sets for which there was insufficient overlap
    locs <- !is.na(score.matrix[,1])
    print(paste("N.gs before overlap prunning:", N.gs))
    N.gs <- sum(locs)
    print(paste("N.gs after overlap prunning:", N.gs))
    if (N.gs > 1) {
        score.matrix <- score.matrix[locs,]
    }
    else {
        # Force vector to matrix if only one gene set specified
        score.matrix <- t(as.matrix(score.matrix[locs,]))
    }
    gs.names <- gs.names[locs]
    gs.descs <- gs.descs[locs]

    # Gene sets with the suffix "_UP" or "_DN" identify gene sets that are known to go up or down across
    # the phenotype of interest. The parameter combine.mode ( = "combine.off", "combine.replace" or "combine.add")
    # controls whether to replace or augment the UP and DN enrichment scores with a single combined score
    # (ES_UP - ES_DN)
    initial.up.entries <- 0
    final.up.entries <- 0
    initial.dn.entries <- 0
    final.dn.entries <- 0
    combined.entries <- 0
    other.entries <- 0

    if (combine.mode == "combine.off") {
        score.matrix.2 <- score.matrix
        gs.names.2 <- gs.names
        gs.descs.2 <- gs.descs
    } else if ((combine.mode == "combine.replace") || (combine.mode == "combine.add")) {
        score.matrix.2 <- NULL
        gs.names.2 <- NULL
        gs.descs.2 <- NULL
        # note: overlap pruning might have reduced N.gs to zero!
        if (N.gs > 0)
        {
            for (i in 1:N.gs) {
                temp <- strsplit(gs.names[i], split="_")
                body <- paste(temp[[1]][seq(1, length(temp[[1]]) -1)], collapse="_")
                suffix <- tail(temp[[1]], 1)
                print(paste("i:", i, "gene set:", gs.names[i], "body:", body, "suffix:", suffix))
                if (suffix == "UP") {  # This is an "UP" gene set
                    initial.up.entries <- initial.up.entries + 1
                    target <- paste(body, "DN", sep="_")
                    loc <- match(target, gs.names)
                    if (!is.na(loc)) { # found corresponding "DN" gene set: create combined entry
                        score <- score.matrix[i,] - score.matrix[loc,]
                        score.matrix.2 <- rbind(score.matrix.2, score)
                        gs.names.2 <- c(gs.names.2, body)
                        gs.descs.2 <- c(gs.descs.2, paste(gs.descs[i], "combined UP & DN"))
                        combined.entries <- combined.entries + 1
                        if (combine.mode == "combine.add") {  # also add the "UP entry
                            score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
                            gs.names.2 <- c(gs.names.2, gs.names[i])
                            gs.descs.2 <- c(gs.descs.2, gs.descs[i])
                            final.up.entries <- final.up.entries + 1
                        }
                    } else { # did not find corresponding "DN" gene set: create "UP" entry
                        score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
                        gs.names.2 <- c(gs.names.2, gs.names[i])
                        gs.descs.2 <- c(gs.descs.2, gs.descs[i])
                        final.up.entries <- final.up.entries + 1
                    }
                } else if (suffix == "DN") { # This is a "DN" gene set
                    initial.dn.entries <- initial.dn.entries + 1
                    target <- paste(body, "UP", sep="_")
                    loc <- match(target, gs.names)
                    if (is.na(loc)) { # did not find corresponding "UP" gene set: create "DN" entry
                        score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
                        gs.names.2 <- c(gs.names.2, gs.names[i])
                        gs.descs.2 <- c(gs.descs.2, gs.descs[i])
                        final.dn.entries <- final.dn.entries + 1
                    } else { # it found corresponding "UP" gene set
                        if (combine.mode == "combine.add") { # create "DN" entry
                            score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
                            gs.names.2 <- c(gs.names.2, gs.names[i])
                            gs.descs.2 <- c(gs.descs.2, gs.descs[i])
                            final.dn.entries <- final.dn.entries + 1
                        }
                    }
                } else { # This is neither "UP nor "DN" gene set: create individual entry
                    score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
                    gs.names.2 <- c(gs.names.2, gs.names[i])
                    gs.descs.2 <- c(gs.descs.2, gs.descs[i])
                    other.entries <- other.entries + 1
                }
            } # end for loop over gene sets
        }

        print(paste("initial.up.entries:", initial.up.entries))
        print(paste("final.up.entries:", final.up.entries))
        print(paste("initial.dn.entries:", initial.dn.entries))
        print(paste("final.dn.entries:", final.dn.entries))
        print(paste("other.entries:", other.entries))
        print(paste("combined.entries:", combined.entries))

        print(paste("total entries:", length(score.matrix.2[,1])))
    }

    if (length(score.matrix.2[,1]) == 0)
    {
        stop("No output gct file written: no gene sets satisfied the min overlap criterion")
    }
    else
    {
        colnames(score.matrix.2) <- sample.names
        row.names(score.matrix.2) <- gs.names.2
        gct <- list(data=score.matrix.2, row.descriptions=gs.descs.2)
        write.gct(gct, output.ds)
    }
} # end of SSGSEA.project.dataset

#
# projects gene expression data onto a single
# gene set by calculating the gene set enrichment score
#
Project.to.GeneSet <- function(
    # data.matrix containing gene expression data
    data.array,

    # exponential weight applied to ranking in calculation of
    # enrichment score
    weight = 0,

    # gene set projecting expression data to
    gene.set)
{

    gene.names <- row.names(data.array)
    n.rows <- dim(data.array)[1]
    n.cols <- dim(data.array)[2]

    ES.vector <- vector(length=n.cols)
    ranked.expression <- vector(length=n.rows, mode="numeric")

    # Compute ES score for signatures in each sample
    for (sample.index in 1:n.cols) {
        # gene.list is permutation (list of row indices) of the normalized expression data, where
        # permutation places expression data in decreasing order
        # Note that in ssGSEA we rank genes by their expression level rather than by a measure of correlation
        # between expression profile and phenotype.
        gene.list <- order(data.array[, sample.index], decreasing=T)

        # gene.set2 contains the indices of the matching genes.
        # Note that when input GCT file is ATARiS-generated, elements of
        # gene.names may not be unique; the following code insures each element
        # of gene.names that is present in the gene.set is referenced in gene.set2
        gene.set2 <- seq(1:length(gene.names))[!is.na(match(gene.names, gene.set))]

        # transform the normalized expression data for a single sample into ranked (in decreasing order)
        # expression values
        if (weight == 0) {
            # don't bother doing calcuation, just set to 1
            ranked.expression <- rep(1, n.rows)
        } else if (weight > 0) {
            # calculate z.score of normalized (e.g., ranked) expression values
            x <- data.array[gene.list, sample.index]
            ranked.expression <- (x - mean(x))/sd(x)
        }

        # tag.indicator flags, within the ranked list of genes, those that are in the gene set
        tag.indicator <- sign(match(gene.list, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag)
        no.tag.indicator <- 1 - tag.indicator
        N <- length(gene.list)
        Nh <- length(gene.set2)
        Nm <-  N - Nh
        # ind are indices into ranked.expression, whose values are in decreasing order, corresonding to
        # genes that are in the gene set
        ind = which(tag.indicator==1)
        ranked.expression <- abs(ranked.expression[ind])^weight

        sum.ranked.expression = sum(ranked.expression)
        # "up" represents the peaks in the mountain plot; i.e., increments in the running-sum
        up = ranked.expression/sum.ranked.expression
        # "gaps" contains the lengths of the gaps between ranked pathway genes
        gaps = (c(ind-1, N) - c(0, ind))
        # "down" contain the valleys in the mountain plot; i.e., the decrements in the running-sum
        down = gaps/Nm
        # calculate the cumulative sums at each of the ranked pathway genes
        RES = cumsum(c(up,up[Nh])-down)
        valleys = RES[1:Nh]-up

        max.ES = max(RES)
        min.ES = min(valleys)

        if( max.ES > -min.ES ){
            arg.ES <- which.max(RES)
        } else{
            arg.ES <- which.min(RES)
        }
        # calculates the area under RES by adding up areas of individual
        # rectangles + triangles
        gaps = gaps+1
        RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES[1:Nh]) - c(valleys,0) ) * (gaps)
        ES = sum(RES)

        ES.vector[sample.index] <- ES

    }
    return(list(ES.vector = ES.vector))

} # end of Project.to.GeneSet


#
# Retrieves a gene set database file (GMT format)
# from GSEA MSigDB. FTP download done by java
# app, which writes db to dest.filename
#
Get.GeneSets.db <- function(
        javaexec,
        jardir,
        gene.sets.db,
        dest.filename)
{
    ssgseaprojection.jar <- paste(jardir, "ssgseaprojection.jar", sep="")
    commons.jar <- paste(jardir, "commons-net-2.0.jar", sep="")

    if (.Platform$OS.type == "windows")
        classpath.sep <- ";"
    else
        classpath.sep <- ":"

    classpath <- paste(".", ssgseaprojection.jar, commons.jar, sep=classpath.sep)

    command <- paste(javaexec, "-cp", classpath ,
                     "org.genepattern.modules.ssgseaprojection.GeneSetsDownloader",
                     gene.sets.db, dest.filename)

    if ((retval <- system(command)) != 0)
        stop("failed to download gene sets db file from GSEA MSigDB; system retval: ", retval)
}

#
# Reads a gene set database file (in GMT file format)
# and creates an R matrix with each row corresponding to a single
# gene set and containing a list of the gene names making up
# that gene set.
#
Read.GeneSets.db <- function(
        gs.db,
        thres.min = 2,
        thres.max = 2000,
        gene.names = NULL)
{

    temp <- readLines(gs.db)
    max.Ng <- length(temp)
    temp.size.G <- vector(length = max.Ng, mode = "numeric")
    for (i in 1:max.Ng) {
        temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
    }
    max.size.G <- max(temp.size.G)
    gs <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
    temp.names <- vector(length = max.Ng, mode = "character")
    temp.desc <- vector(length = max.Ng, mode = "character")
    gs.count <- 1
    for (i in 1:max.Ng) {
        gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
        gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
        gene.set.name <- gs.line[1]
        gene.set.desc <- gs.line[2]
        gene.set.tags <- vector(length = gene.set.size, mode = "character")
        for (j in 1:gene.set.size) {
            gene.set.tags[j] <- gs.line[j + 2]
        }
        if (is.null(gene.names)) {
            existing.set <- rep(TRUE, length(gene.set.tags))
        } else {
            existing.set <- is.element(gene.set.tags, gene.names)
        }
        set.size <- length(existing.set[existing.set == T])
        if ((set.size < thres.min) || (set.size > thres.max)) next
        temp.size.G[gs.count] <- set.size
        gs[gs.count,] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
        temp.names[gs.count] <- gene.set.name
        temp.desc[gs.count] <- gene.set.desc
        gs.count <- gs.count + 1
    }
    Ng <- gs.count - 1
    gs.names <- vector(length = Ng, mode = "character")
    gs.desc <- vector(length = Ng, mode = "character")
    size.G <- vector(length = Ng, mode = "numeric")

    gs.names <- temp.names[1:Ng]
    gs.desc <- temp.desc[1:Ng]
    size.G <- temp.size.G[1:Ng]

    return(list(N.gs = Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc, size.G = size.G, max.N.gs = max.Ng))
} # end of Read.GeneSets.db



