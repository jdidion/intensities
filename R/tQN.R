# Perform tQN normalization of a set of samples. Requires a cluster file generated from
# reference samples, which can be created with the create.cluster.file() function.
# Requires: bioconductor
# This script is adapted from code provided by Johan Staaf.
# Staaf et al., BMC Bioinformatics, 2008, doi:10.1186/1471-2105-9-409.

library(multicore)

tryCatch(library(limma), error=function(e) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("limma")
})

is.dir <- function(d) file.exists(d) && file.info(d)$isdir

is.num <- function(x) !(is.na(x) | is.nan(x) | is.infinite(x))

theta <- function(data) 2/pi*atan(data$Y/data$X)

bound <- function(x, low, high) {
    n <- is.num(x)
    x[n & x < low] <- low
    x[n & x > high] <- high
    x
}

my.mean <- function(x, minlen) {
    x <- x[is.num(x)]
    l <- length(x)
    if (l < minlen) {
        NA
    }
    else if (l == 1) {
        x
    }
    else {
        mean(x)
    }
}

my.sd <- function(x, minlen) {
    x <- x[is.num(x)]
    l <- length(x)
    if (l < minlen) {
        NA
    }
    else if (l == 1) {
        0
    }
    else {
        sd(x)
    }
}

# Prepare a data frame with reporter annotations from a file
#
# reps: a tab-delim file or data frame with at least 5 columns: reporterId, Chr, Position, X.allele and Y.allele,
# where alleles are specified as A/C/G/T.
prepare.reporters <- function(reps) {
    if (is.character(reps)) {
        reps <- my.load(reps)
    }
    if (ncol(reps) < 5) {
        stop("Invalid reps table")
    }
    rownames(reps) <- reps$reporterId
    # code sex chromosomes and mitochondria as numeric
    reps$Chr[reps$Chr == "X"] <- 23
    reps$Chr[reps$Chr == "Y"] <- 24
    reps$Chr[reps$Chr == "M"] <- 25
    reps$Chr <- as.numeric(reps$Chr)
    # Filter unannotated reporters
    reps <- reps[!is.na(reps$Chr) & reps$Chr > 0 & !is.na(reps$Position) & reps$Position > 0,]
    # Sort reporters by chromosome and position
    reps[order(reps$Chr, reps$Position),]
}

# Normalize reference sample intensities prior to creating a cluster file.
#
# reps: a data frame of reporters (loaded/prepared using prepare.reporters).
# data.file: path of a tab-delimited data file that must have the following columns (including a header with 
# these exact names): reporterId, X, Y and EITHER call OR Allele1.Forward and Allele2.Forward. The formats are: 
# allele1 (A/B/H/N), allele2 (AA/BB/AB/NC), nucleotide1 (A/C/G/T/H/N), or nucleotide2 (AA/CC/GG/TT/AC/AG.../NC). 
# By default, call is assumed to be in allele2 format and the allele columns are assumed to be nucleotides (so 
# the default format would be nucleotide2).
# QN.thresholds: thresholds for the X and Y values after quantile normalization.
# save: whether to save the data frame to a tab-delimited text file
# output.dir: ignored if save is FALSE. Otherwise, if NULL, the original files will be overwritten, otherwise 
# the file will be written with the same name but to the specified directory. The new files contain the same data 
# but with extra columns, some rows may be removed, and the genotype call format may be changed.
tqn.normalize <- function(reps, data.file, call.format="allele1", QN.thresholds=c(1.5,1.5), save=TRUE, output.dir=NULL) {   
    reps <- prepare.reporters(reps)

    data <- read.table(data.file, sep="\t", stringsAsFactors=FALSE, header=TRUE, na.strings=c("", "NA"))
    rownames(data) <- data$reporterId
    
    # remove unused genotypes
    data <- data[rownames(reps),]

    # convert to allele2 format
    if ("call" %in% colnames(data)) {
        calls <- data$call
        call.format = match.arg(call.format, c("allele2", "allele1", "nucleotide1", "nucleotide2"))
        if (call.format == "nucleotide2") {
            nucs <- do.call(rbind, strsplit(calls, ''))
        }
    }
    else {
        call.format = match.arg(call.format, c("nucleotide2", "allele1", "allele2", "nucleotide1"))
        if (call.format == "allele2") {
            calls <- paste(data$Allele1.Forward, data$Allele2.Forward, sep="")
        }
        else {
            nucs <- data[,c("Allele1.Forward","Allele2.Forward")]
        }           
    }

    if (call.format == "allele1") {
        calls[calls == "A"] <- "AA"
        calls[calls == "B"] <- "BB"
        calls[calls == "H"] <- "AB"
        calls[calls == "N"] <- "NC"
    }
    else if (call.format == "nucleotide1") {
        AA <- calls == reps$X.allele
        BB <- calls == reps$Y.allele
        NC <- calls == "N"
        calls[AA] <- "AA"
        calls[BB] <- "BB"
        calls[NC] <- "NC"
        calls[!(AA|BB|NC)] <- "AB"
    }
    else if (call.format == "nucleotide2") {
        NC <- nucs[,1] == "N"
        hom <- nucs[,1] == nucs[,2]
        calls[!NC & hom & nucs[,1] == "X.allele"] <- "AA"
        calls[!NC & hom & nucs[,1] == "Y.allele"] <- "BB"
        calls[!NC & !hom] <- "AB"
        calls[NC] <- "NC"
    }

    data$call = calls

    # quantile normalize intensities
    inorm <- normalizeQuantiles(data[,c("X","Y")])
    colnames(inorm) <- c("X", "Y")
    
    # Normalized X, Y and R (without thresholding)
    data$XCorrected <- inorm$X
    data$YCorrected <- inorm$Y
    data$RCorrected <- inorm$X + inorm$Y

    na <- data$X <= 0 | data$Y <= 0

    # Thesholding of QN effect
    aff.x <- !na & (inorm$X / data$X) > QN.thresholds[1]
    inorm$X[aff.x] <- QN.thresholds[1] * data$X[aff.x]

    aff.y <- !na & (inorm$Y / data$Y) > QN.thresholds[2]
    inorm$Y[aff.y] <- QN.thresholds[2] * data$Y[aff.y]

    # Corrected T
    T <- theta(inorm)
    T[na] <- NA
    T[!na] <- bound(T[!na], 0 ,1)
    data$TCorrected <- T

    if (save) {
        if (is.null(output.dir)) {
            outfile <- data.file
        }
        else {
            if (!file.exists(output.dir)) {
                dir.create(output.dir)
            }
            outfile <- file.path(output.dir, basename(data.file))
        }
        write.table(data, outfile, sep="\t", quote=FALSE, append=FALSE, row.names=FALSE, col.names=TRUE)
    }
    
    invisible(data)
}

# Create a tQN cluster file from a set of reference samples. 
#
# The memory requirement is ~ M*N*6*B, where M is the number of reporters, N is the number of samples and B is 
# the size of an integer, in bytes (4 for 32-bit R and 8 for 64-bit R). Therefore, if there are a large number 
# of reporters and/or samples, you will need to first run the normaliztion over all the files, then run this 
# function using a chunk size that corresponds to your available memory.
#
# The files argument is either a vector of file names, a directory name, or the name of a file that lists all of
# the data files. This function expects that the last five columns in each data file will be: call, XCorrected, 
# YCorrected, RCorrected, TCorrected (as this is the output of the tqn.normalize function). If this isn't true, 
# you must supply the col.idxs parameter, which is a vector of indexes for the call, RCorrected and TCorrected columns.
#
# The extra arguments are passed to tqn.normalize.
create.cluster.file <- function(reps, ref.files, cluster.file, pre.normalized=FALSE, chunk.size=NULL, col.idxs=NULL, 
        min.samples.per.cluster=1, ...) {
    reps <- prepare.reporters(reps)

    if (length(ref.files) == 1) {
        if (is.dir(ref.files[1])) {
            ref.files <- sapply(list.files(ref.files), function(f) file.path(ref.files, f))
        }
        else {
            ref.files <- scan(ref.files, "character")
        }
    }
    nfiles <- length(ref.files)

    split <- !is.null(chunk.size)
    if (split && !pre.normalized) {
        stop("Processing files in chunks requires data to be pre-normalized")
    }

    chunk <- 0
    done <- FALSE
    cols <- c("AA_T","AB_T","BB_T","AA_R","AB_R","BB_R")

    while (!done) {
        if (split) {
            r <- reps[((chunk * chunk.size) + 1):min((chunk + 1) * chunk.size, nrow(reps)),]
            if (nrow(r) < chunk.size) {
                done <- TRUE
            }
        }
        else {
            r <- reps
            done <- TRUE
        }

        # Aggregate R and T values by allele using a 3D array.
        pop.data <- array(NA, dim=c(nrow(r), nfiles, 6), dimnames=list(rownames(r), 1:nfiles, cols))
        
        for (i in 1:nfiles) {
            if (pre.normalized) {
                data <- read.table(ref.files[i], sep="\t", stringsAsFactors=FALSE, na.strings=c("", "NA", "NaN"), header=chunk == 0, 
                    skip=ifelse(split, 0, chunk * chunk.size), nrows=ifelse(split, -1, chunk.size))
            }
            else {
                data <- tqn.normalize(reps, ref.files[i], ...)
            }
            
            ixs <- col.idxs
            if (is.null(ixs)) {
                n <- ncol(data)
                ixs <- c(n-4, n-1, n)
            }

            AA <- data[,ixs[1]] == "AA"
            pop.data[AA, i, c("AA_R", "AA_T")] <- as.matrix(data[AA, ixs[2:3]])
            
            AB <- data[,ixs[1]] == "AB"
            pop.data[AB, i, c("AB_R", "AB_T")] <- as.matrix(data[AB, ixs[2:3]])

            BB <- data[,ixs[1]] == "BB"
            pop.data[BB, i, c("BB_R", "BB_T")] <- as.matrix(data[BB, ixs[2:3]])
        }

        means <- t(apply(pop.data, 1, apply, 2, my.mean, min.samples.per.cluster))
        devs <- t(apply(pop.data, 1, apply, 2, my.sd, min.samples.per.cluster))

        write.table(data.frame(reporterId=rownames(reps), do.call(cbind, lapply(cols, function(x) {
                tab <- cbind(means[,x], devs[,x])
                colnames(tab) <- c(paste(x, "Mean", sep="_"), paste(x, "Dev", sep="_"))
                tab
            })), stringsAsFactors=FALSE), 
            cluster.file, sep="\t", append=split && chunk > 0, quote=FALSE, row.names=FALSE, col.names=TRUE)

        chunk <- chunk + 1
    }

    nfiles
}       

prepare.clusters <- function(clusters) {
    if (is.character(clusters)) {
        clusters <- read.table(clusters, sep="\t", stringsAsFactors=FALSE, header=TRUE, na.strings=c("","NA","NaN"))
    }
    rownames(clusters) <- clusters$reporterId
    clusters
}

# Run tQN normalization for a set of samples. sample.files may be a directory or a vector of files. If 
# output.dir is defined, updated files will be written there, otherwise the original files will be 
# overwritten. Output is in native tQN format. A perl script is provided with the tQN package that can 
# translate this into other formats.
tQN.files <- function(reps, clusters, sample.files, output.dir=NULL, QN.thresholds=c(1.5,1.5), mc.cores=14) {
    reps <- prepare.reporters(reps)
    clusters <- prepare.clusters(clusters)

    if (length(sample.files) == 1 && is.dir(sample.files[1])) {
        sample.files <- sapply(list.files(sample.files), function(f) file.path(sample.files, f))
    }

    mclapply(sample.files, function(f) {
        data <- read.table(f, sep="\t", stringsAsFactors=FALSE, header=TRUE, na.strings=c("","NA","NaN"))
        idcol <- ifelse("reporterId" %in% colnames(data), "reporterId", "Name")
        rownames(data) <- data[,idcol]

        norm.data <- tQN(data, reps, clusters, QN.thresholds)
        
        # Write output
        if (is.null(output.dir)) {
            outfile <- f
        }
        else {
            if (!file.exists(output.dir)) {
                dir.create(output.dir)
            }
            outfile <- file.path(output.dir, basename(f))
        }

        write.table(norm.data,
            outfile, sep="\t", quote=FALSE, append=FALSE, row.names=FALSE,
            col.names=c("Name", "Chr", "Position", "tQN B Allele Frequency", "tQN Log R Ratio", "tQN X", "tQN Y"))

        outfile
    }, mc.cores=mc.cores)
}

tQN.database <- function(reps="~/data/arrays/megamuga/tQN/mm_reporters.RData",
        clusters="~/data/arrays/megamuga/tQN/mm_tQN_clusters.txt", sampleIDs, 
        database="~/data/arrays/megamuga/megamuga.db", ...) {
    library(RSQLite)
    db <- dbConnect(dbDriver("SQLite"), dbname=database)
    mat <- dbGetPreparedQuery(db, "SELECT sampleID, snpID, X, Y FROM genotypes WHERE sampleID=?", data.frame(sampleID=sampleIDs))
    dbDisconnect(db)
    tQN.matrix(reps, clusters, mat, ...)
}

# Run tQN normalization for a set of samples in a matrix (sampleID, snpID, X, Y). If output.dir is
# not null, result files will be written there as <sampleID>.norm.txt.
tQN.matrix <- function(reps, clusters, mat, output.dir=NULL, QN.thresholds=c(1.5,1.5), mc.cores=14) {
    reps <- prepare.reporters(reps)
    clusters <- prepare.clusters(clusters)

    ids <- unique(mat[,1])
    mclapply(ids, function(i) {
        w <- mat[,1] == i
        data <- mat[w,2:4]
        rownames(data) <- data[,1]
        colnames(data) <- c("reporterId","X","Y")
        
        norm.data <- tQN(data, reps, clusters, QN.thresholds)

        # Write output
        if (is.null(output.dir)) {
            norm.data
        }
        else {
            if (!file.exists(output.dir)) {
                dir.create(output.dir)
            }
            outfile <- file.path(output.dir, paste(i, ".norm.txt", sep=""))

            write.table(norm.data,
                outfile, sep="\t", quote=FALSE, append=FALSE, row.names=FALSE,
                col.names=c("Name", "Chr", "Position", "tQN B Allele Frequency", "tQN Log R Ratio", "tQN X", "tQN Y"))

            outfile
        }
    }, mc.cores=mc.cores)
}

tQN <- function(data, reps, clusters, QN.thresholds=c(1.5,1.5)) {
    # Remove any rows that aren't in the cluster file
    i <- intersect(rownames(clusters), rownames(data))
    data <- data[i,]
    ref <- clusters[i,]

    # Quantile normalization
    inorm <- normalizeQuantiles(data[,c("X","Y")])

    # Calculate R
    R <- inorm$X + inorm$Y

    na <- data$X <= 0 | data$Y <= 0

    # Thesholding of QN effect
    aff.x <- !na & (inorm$X / data$X) > QN.thresholds[1]
    inorm$X[aff.x] <- QN.thresholds[1] * data$X[aff.x]

    aff.y <- !na & (inorm$Y / data$Y) > QN.thresholds[2]
    inorm$Y[aff.y] <- QN.thresholds[2] * data$Y[aff.y]

    T <- theta(inorm)
    T[na] <- NA
    T[!na] <- bound(T[!na], 0 ,1)
    
    # Calculate tQN X and Y to fit theta and R
    Y <- R * tan(T * pi/2) / (1 + tan(T * pi/2))
    X <- R - Y

    # Calculate BAF and LRR from corrected R and T
    baflrr <- baf.lrr(R, T, ref)

    data.frame(i, reps[i, c("Chr", "Position")], baflrr, X, Y, stringsAsFactors=FALSE)
}

# Compute BAF
baf.lrr <- function(R, T, ref) {
    BAF <- T
    LRR <- R
    
    med.tAA <- median(ref$AA_T_Mean, na.rm=TRUE)
    med.tBB <- median(ref$BB_T_Mean, na.rm=TRUE)
    
    for (i in 1:length(T)) {
        if (!is.num(T[i])) {
            BAF[i] <- NaN
            LRR[i] <- NaN
            next
        }
    
        th <- T[i]
        rr <- R[i]
        rAA <- ref$AA_R_Mean[i]
        rAB <- ref$AB_R_Mean[i]
        rBB <- ref$BB_R_Mean[i]
        tAA <- ref$AA_T_Mean[i]
        tAB <- ref$AB_T_Mean[i]
        tBB <- ref$BB_T_Mean[i]

        e.tAA <- is.num(tAA)
        e.tAB <- is.num(tAB)
        e.tBB <- is.num(tBB)
    
        # 0: Test for inconsistencies between tAA/tAB/tBB and rAA/rAB/rBB
        if (((e.tAA & e.tAB) && tAA > tAB) ||
                ((e.tAA & e.tBB) && tAA > tBB) ||
                ((e.tAB & e.tBB) && tAB > tBB)) {
            BAF[i] <- LRR[i] <- NaN
        }
        # 1: Triple blank SNP
        else if (!(e.tAA|e.tAB|e.tBB)) {
            BAF[i] <- LRR[i] <- NaN
        }
        # 2: Blank for AB, AA, while positive for BB
        else if (!(e.tAA|e.tAB) & e.tBB) {
            if (th >= tBB) {
                BAF[i] <- 1
                LRR[i] <- ifelse(rBB[i] <= 0, NaN, log2(rr / rBB[i]))
            }
            else {
                BAF[i] <- LRR[i] <- NaN
            }
        }
        # 3: Blank for AB, BB, while positive for AA
        else if (e.tAA & !(e.tAB|e.tBB)) {
            if (th <= tAA) {
                BAF[i] <- 0
                LRR[i] <- ifelse(tAA[i] <= 0, NaN, log2(rr / tAA[i]))
            }
            else {
                BAF[i] <- LRR[i] <- NaN
            }
        }
        # 4: Blank for AB while positive for AA & BB
        else if (e.tAA & !e.tAB & e.tBB) {
            # no AB cluster exist for this SNP, while AA & BB exists. Set it to the closest of AA or BB
            min.index <- which.min(c(abs(tAA - th), abs(tBB - th)))
            if (min.index == 1 && th < tAA) {
                BAF[i] <- 0
                LRR[i] <- ifelse(tAA[i] <= 0, NaN, log2(rr / tAA[i]))
            }
            else if (min.index != 1 && th >= tBB) {
                BAF[i] <- 1
                LRR[i] <- ifelse(rBB[i] <= 0, NaN, log2(rr / rBB[i]))
            }
            else {
                BAF[i] <- LRR[i] <- NaN
            }
        }
        # 5: Blank for AA while positive for AB & BB
        else if (!e.tAA & e.tAB & e.tBB) {
            if (th >= tBB) {
                BAF[i] <- 1
                LRR[i] <- NaN
            }
            # 5.1: SNP is "correctly between" ref$AB_T_Mean and ref$BB_T_Mean
            else if (th >= tAB) {
                # interpolate as SNP is expected to be between ref$AB_T_Mean and ref$BB_T_Mean
                BAF[i] <- 0.5 + 0.5 * (th - tAB) / (tBB - tAB)
                eR <- rAB + ((th - tAB) * (rBB - rAB) / (tBB - tAB))
                LRR[i] <- ifelse(eR <= 0, NaN, log2(rr / eR))
            }
            # 5.2: Heterozygous SNP is subjected to deletion or UPD of allele B making it unexectedly to be 
            # between ref$AA_T_Mean and ref$AB_T_Mean where it normally should not NOT BE.
            else {
                BAF[i] <- ifelse(th < med.tAA, 0, 0.5 * (th - med.tAA) / (tAB - med.tAA))
                LRR[i] <- NaN                   
            }
        }
        # 6: Blank for BB while positive for AA & AB
        else if (e.tAA & e.tAB & !e.tBB) {
            if (th < tAA) {
                BAF[i] <- 0
                LRR[i] <- NaN
            }
            # 6.1: SNP is "correctly between" ref$AA_T_Mean and ref$AB_T_Mean
            else if (th <= tAB) {
                #interpolate as SNP is expected to be between ref$AB_T_Mean and ref$BB_T_Mean
                BAF[i] <- 0.5* (th - tAA) / (tAB - tAA)
                eR <- rAA + ((th - tAA) * (rAB - rAA) / (tAB - tAA))
                LRR[i] <- ifelse(eR <= 0, NaN, log2(rr / eR))
            }
            # 2: Heterozygous SNP is subjected to deletion or UPD of allele A making it unexectedly to be 
            # between ref$AB_T_Mean and ref$BB_T_Mean where it normally should not NOT BE.
            else {
                BAF[i] <- ifelse(th > med.tBB, 1, 0.5 + 0.5 * (th - tAB) / (med.tBB[i] - tAB))
                LRR[i] <- NaN
            }
        }
        # 7: positive for AA & BB & AB, Illumina style calculation
        else {
            if (th < tAB) {
                BAF[i] <- ifelse(th < tAA, 0, 0.5 * (th - tAA) / (tAB - tAA))
                eR <- rAA + ((th - tAA) * (rAB - rAA) / (tAB - tAA))
            }
            else {
                BAF[i] <- ifelse(th >= tBB, 1, 0.5 + 0.5 * (th - tAB) / (tBB - tAB))
                eR <- rAB + ((th - tAB) * (rBB - rAB) / (tBB - tAB))
            }
            LRR[i] <- ifelse(eR <= 0, NaN, log2(rr / eR))
        }
    }

    cbind(BAF=bound(BAF, 0, 1), LRR=LRR)
}

