# © 2012 Broad Institute, Inc.  All rights reserved.

# ssGSEA Projection

# processes the cmd line for the ssGSEA.project.dataset
ssGSEA.projection.cmdline <- function(...)
{
    input.gct.filename               <- NA
    output.prefix                    <- NA
    gene.sets.database               <- NA
    gene.sets.db.filename            <- NA
    gene.sets.db.list.filename       <- NA
    gene.symbol.column.name          <- "Name"
    gene.set.selection               <- "ALL"
    sample.normalization.method      <- "rank"
    weighting.exponent               <- 0.75
    min.overlap                      <- 1

    args <- list(...)
    for (i in 1:length(args))
    {
        flag <- substring(args[[i]], 1, 2)
        value <- substring(args[[i]], 3, nchar(args[i]))

        if (value == '')
        {
            next
        }
        else if (flag == '-l')
        {
            libdir <- value
        }
        else if (flag == '-j')
        {
            javaexec <- value
        }
        else if (flag == '-i')
        {
            input.gct.filename <- value
        }
        else if (flag == '-o')
        {
            output.prefix <- value
        }
        else if (flag == '-G')
        {
            gene.sets.database <- value
        }
        else if (flag == '-d')
        {
            gene.sets.db.filename <- value
        }
        else if (flag == '-D')
        {
            gene.sets.db.list.filename <- value
        }
        else if (flag == '-c')
        {
            gene.symbol.column.name <- value
        }
        else if (flag == '-s')
        {
            gene.set.selection <- unlist(strsplit(value,','))
        }
        else if (flag == '-n')
        {
            sample.normalization.method <- value
        }
        else if (flag == '-w')
        {
            weighting.exponent <- as.numeric(value)
        }
        else if (flag == '-v')
        {
            min.overlap <- as.integer(value)
        }
        else
            stop("Unknown option", flag)

    }

    if (is.na(input.gct.filename))
        stop("Missing input.gct.filename")

    if (is.na(output.prefix))
    {
        temp <- strsplit(input.gct.filename, split="/") # Extract input file name
        s <- length(temp[[1]])
        input.file.name <- temp[[1]][s]
        temp <- strsplit(input.file.name, split=".gct")
        output.prefix <-  paste(temp[[1]][1],".PROJ", sep="")
    }

    gene.sets.dbfile.list <- NA
    if (!is.na(gene.sets.db.filename))
        gene.sets.dbfile.list <- c(gene.sets.db.filename)

    if (!is.na(gene.sets.db.list.filename))
        gene.sets.dbfile.list <- readLines(gene.sets.db.list.filename)

    if (length(gene.sets.dbfile.list) == 0)
        gene.sets.dbfile.list <- NA

    setup(libdir)
    source(paste(libdir,"ssGSEAProjection.Library.R", sep=''))

    suppressWarnings(ssGSEA.project.dataset(javaexec,
                                            libdir,
                                            input.gct.filename,
                                            paste(output.prefix, ".gct", sep=""),
                                            gene.sets.database,
                                            gene.sets.dbfile.list = gene.sets.dbfile.list,
                                            gene.symbol.column = gene.symbol.column.name,
                                            gene.set.selection = gene.set.selection,
                                            sample.norm.type = sample.normalization.method,
                                            weight = weighting.exponent,
                                            min.overlap = min.overlap))
}

setup <- function(libdir)
{
    source(paste(libdir, "common.R", sep=''))
    setLibPath(libdir)
    install.required.packages(libdir)
}

install.required.packages <- function(libdir)
{
    info(libdir)
    # no non-base packages required by this module
}


