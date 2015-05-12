# Methods for detecting copy expansion or deletion, relative to a known karyotypically normal
# sample. Metric calculations are adapted from Peiffer et al., Genome Research, 2006,
# doi:10.1101/gr.5402306. Unlike the validation assay, these methods require genotype data in
# a specific format (A, B, H or N calls).

# High-level function that takes data in standard format and calculates paramters from canonical
# samples (if any) and then calculates copy number metrics for test samples (if any). If the params
# argument is specified, any canonical samples are ignored. Each input parameter can be a matrix or
# delmited file of the specified format:
# samples: sampleID, name, sex, and whether or not the sample is canonical (TRUE/FALSE), i.e.
# whether it should be used to calculate the genotype cluster parameters. If the fourth column
# is not specified, all samples are assumed to be canonical if params is NULL, otherwise all
# samples are assumed to be non-canonical.
# snps: snpID, chrm, pos. Should be sorted by chromosome and position.
# data: sampleID, snpID, x, y, call where x and y are the normalized intensity values. If there
# are only four columns, then genotype calls will be taken from the genotypes matrix.
# genotypes: Rows are SNPs and columns are samples.
# params: a params list previously obtained from calc.params().
get.copy.number.metrics <- function(samples, snps, data, genotypes=NULL, params=NULL, sep="\t",
        plot=FALSE, ...) {
    # Prepare sample data
    print("Loading samples")
    if (is.character(samples)) {
        samples <- read.table(samples, header=TRUE, sep=sep, row.names=1, stringsAsFactors=FALSE,
            comment.char="", quote="")
    }
    else if (!is.data.frame(samples)) {
        samples <- as.data.frame(samples, stringsAsFactors=FALSE)
    }
    colnames(samples)[1:2] <- c("name", "sex")
    if (ncol(samples) == 2) {
        samples$canonical = is.null(params)
    }
    else {
        colnames(samples)[3] <- "canonical"
    }

    # Prepare SNP data
    print("Loading SNPs")
    if (is.character(snps)) {
        snps <- read.table(snps, header=TRUE, sep=sep, row.names=1, stringsAsFactors=FALSE,
            comment.char="", quote="")
    }
    else if (!is.data.frame(snps)) {
        snps <- as.data.frame(snps, stringsAsFactors=FALSE)
    }
    colnames(snps) <- c("chrm", "pos")
    snps$chrm <- as.character(snps.chrm)

    # Prepare intensity/genotype data
    print("Loading intensity")
    if (is.character(data)) {
        data <- read.table(data, header=TRUE, sep=sep, row.names=NULL, stringsAsFactors=FALSE,
            comment.char="", quote="",
            colClasses=c("character", "character", "numeric", "numeric", "character"))
    }
    else if (!is.data.frame(data)) {
        data <- as.data.frame(data, stringsAsFactors=FALSE)
    }
    colnames(data)[1:4] <- c("sampleID", "snpID", "x", "y")
    if (ncol(data) == 4 && !is.null(genotypes)) {
        if (is.character(genotypes)) {
            genotypes <- read.table(genotypes, header=TRUE, sep=sep, row.names=1, stringsAsFactors=FALSE)
        }

        require(reshape2)
        genotypes <- melt(genotypes)
        colnames(genotypes) <- c("snpID", "sampleID", "call")

        m <- match(paste(data$sampleID, data$snpID), paste(genotypes$sampleID, genotypes$snpID))
        data$call <- genotypes[m, "call"]
    }
    else {
        colnames(data)[5] <- "call"
    }

    # Calculate parameters from canonical samples if necessary
    if (is.null(params)) {
        print("Calculating reference parameters")
        canonical.samples <- samples[samples$canonical,]
        # Convert SNPs to list format
        params <- calc.params(canonical.samples, as.snp.list(snps),
            data[!is.na(match(data$sampleID, canonical.samples$sampleID)),])
    }

    if (is.null(params)) {
        stop("Missing required parameters from canonical samples")
    }

    result <- list(samples=samples, snps=snps, params=params)

    # Calcluate metrics for test samples
    test <- !samples$canonical
    if (any(test)) {
        print("Calculating metrics")
        test.samples <- samples[test,]
        geno.list <- as.geno.list(data[!is.na(match(data$sampleID, test.samples$sampleID)),])
        result$metrics = calc.copy.number.metrics(test.samples, params, geno.list)

        if (plot) {
            plot.copy.number.metrics(test.samples, snps, result$metrics, ...)
        }
    }

    result
}

## Low-level functions. These use the following data types:
# snps.list: A list of three items: autosomal, X, and Y chromosome snps. Each item is a two-column
# matrix of chrm, pos, with snpIDs as row names.
# samples: A two-column matrix: name, sex (m/f), with sample IDs as rownames.
# geno: A five-column matrix: sampleID, snpID, x, y, call, where call is A/B/H/N.
# geno.list: A list of geno matricies, one per sample. Each matrix only needs three columns:
# call, x, y, with snpIDs as rownames.
# params: A list of four items: autosomal, male X, female X and Y chromosome parameters. Each item
# is a 13-column matrix: thetaAA, thetaAB, thetaBB, rAA, rAB, rBB, sAA, sAB, sBB, mrA, brA, mrB, brB,
# msA, bsA, mrB, brB,
# where AA, AB and BB are the centroids of "canonical" clusters, r and s are the mean and standard
# deviation of R for each cluster, and the m and b values are the slopes and intercepts of lines
# that pass through (AA,AB) and (AB,BB). Row names are snpIDs.
# metrics: A three-column metric matrix for: call, baf and lrr, with snpIDs as row names, and where
# baf is the B Allele Frequency and lrr is the Log R Ratio (both defined in Peiffer et al).
# metrics.list: A list of metric matricies, one per sample.

# Convert a single SNP matrix into list format.
as.snp.list <- function(snps) {
    X <- snps$chrm == "X"
    Y <- snps$chrm == "Y"
    auto <- !(X|Y)
    list(auto=snps[auto,], X=snps[X,], Y=snps[Y,])
}

# Convert a single geno matrix to list format.
as.geno.list <- function(geno) {
    by(geno, geno$sampleID, function(g) data.frame(g[,c("x","y","call")], row.names=g$snpID))
}

# Calculate control parameters based on intensity (normalized x and y) values from
# a set of canonical samples.
calc.params <- function(samples, snps.list, geno, min.samples.per.cluster=5) {
    m.ids <- rownames(samples)[samples$sex == 'm']
    m <- geno[!is.na(match(geno$sampleID, m.ids)),]
    f.ids <- rownames(samples)[samples$sex == 'f']
    f <- geno[!is.na(match(geno$sampleID, f.ids)),]
    list(auto=.calc.params(snps.list$auto, geno, min.samples.per.cluster, TRUE),
        mX=.calc.params(snps.list$X, m, min.samples.per.cluster, FALSE),
        fX=.calc.params(snps.list$X, f, min.samples.per.cluster, TRUE),
        Y=.calc.params(snps.list$Y, m, min.samples.per.cluster, FALSE))

}

.calc.params <- function(snps, geno, min.samples.per.cluster, require.hets=TRUE) {
    if (is.null(snps) || nrow(snps) == 0) return(NULL)
    
    ## Filter out genotypes/SNPs that don't meet criteria
    # Remove genotypes outside the valid range of 0 <= x,y <= 1
    geno <- geno[geno$x >= 0 & geno$y >= 0,]
    # Remove SNPs with too few samples in each cluster
    counts <- as.data.frame(table(geno$snpID, geno$call)[rownames(snps), c("A","H","B")])
    valid <- (counts$A >= min.samples.per.cluster) & (counts$B >= min.samples.per.cluster)
    if (require.hets) {
        valid <- valid & (counts$H >= min.samples.per.cluster)
    }
    snps <- snps[valid,]
    if (sum(valid) == 0) {
        return(NULL)
    }

    if (!require.hets) {
        geno <- geno[geno$call != "H",]
    }

    i <- intersect(rownames(snps), geno$snpID)
    snps <- snps[i,,drop=FALSE]
    geno <- geno[!is.na(match(geno$snpID, i)),]

    theta <- (2 / pi) * atan(geno$y / geno$x)
    theta[is.nan(theta)] <- 0
    r <- geno$x + geno$y

    # Calculate the mean theta for each genotype and each SNP. If require.hets is TRUE,
    # we are only interested in SNPs that have all three genotype clusters, otherwise
    # they must have at least A and B clusters.
    t.mean <- as.data.frame(force.matrix(tapply(theta,
        list(geno$snpID, factor(geno$call, levels=c("A","H","B","N"))), mean)[rownames(snps), 1:3],
        default.name=i))
    colnames(t.mean) <- c("thetaAA", "thetaAB", "thetaBB")

    # Calculate the mean R for each genotype and each SNP and reduce to only
    # the SNPs with valid theta values.
    r.mean <- as.data.frame(force.matrix(tapply(r,
        list(geno$snpID, factor(geno$call, levels=c("A","H","B","N"))), mean)[rownames(snps), 1:3],
        default.name=i))
    colnames(r.mean) <- c("rAA", "rAB", "rBB")
    r.sd <- as.data.frame(force.matrix(tapply(r,
        list(geno$snpID, factor(geno$call, levels=c("A","H","B","N"))), sd)[rownames(snps), 1:3],
        default.name=i))
    colnames(r.sd) <- c("sAA", "sAB", "sBB")

    # Calculate slope (m) and intercept (b) for the mean (r) and standard deviation (s) for each
    # SNP for (AA,AB) and (AB,BB). If require.hets is FALSE, then m and b will only be calculated
    # for (AA,BB) and mrB,brB will be NA.
    params <- data.frame(t.mean, r.mean, r.sd, mrA=NA, brA=NA, mrB=NA, brB=NA, msA=NA, bsA=NA, msB=NA, bsB=NA)
    het.missing <- is.na(params$thetaAB)
    if (any(het.missing)) {
        params[het.missing, "mrA"] <- (params[het.missing, "rBB"] - params[het.missing, "rAA"]) / (params[het.missing, "thetaBB"] - params[het.missing, "thetaAA"])
        params[het.missing, "brA"] <- params[het.missing, "rAA"] - (params[het.missing, "mrA"] * params[het.missing, "thetaAA"])
        params[het.missing, "msA"] <- (params[het.missing, "sBB"] - params[het.missing, "sAA"]) / (params[het.missing, "thetaBB"] - params[het.missing, "thetaAA"])
        params[het.missing, "bsA"] <- params[het.missing, "sAA"] - (params[het.missing, "msA"] * params[het.missing, "thetaAA"])
    }
    params[!het.missing, "mrA"] <- (params[!het.missing, "rAB"] - params[!het.missing, "rAA"]) / (params[!het.missing, "thetaAB"] - params[!het.missing, "thetaAA"])
    params[!het.missing, "brA"] <- params[!het.missing, "rAA"] - (params[!het.missing, "mrA"] * params[!het.missing, "thetaAA"])
    params[!het.missing, "msA"] <- (params[!het.missing, "sAB"] - params[!het.missing, "sAA"]) / (params[!het.missing, "thetaAB"] - params[!het.missing, "thetaAA"])
    params[!het.missing, "bsA"] <- params[!het.missing, "sAA"] - (params[!het.missing, "msA"] * params[!het.missing, "thetaAA"])
    params[!het.missing, "mrB"] <- (params[!het.missing, "rBB"] - params[!het.missing, "rAB"]) / (params[!het.missing, "thetaBB"] - params[!het.missing, "thetaAB"])
    params[!het.missing, "brB"] <- params[!het.missing, "rAB"] - (params[!het.missing, "mrB"] * params[!het.missing, "thetaAB"])
    params[!het.missing, "msB"] <- (params[!het.missing, "sBB"] - params[!het.missing, "sAB"]) / (params[!het.missing, "thetaBB"] - params[!het.missing, "thetaAB"])
    params[!het.missing, "bsB"] <- params[!het.missing, "sAB"] - (params[!het.missing, "msB"] * params[!het.missing, "thetaAB"])

    params
}

# Calculate the metrics for copy number analysis. SNPs with noisy intensity data
# can be eliminated using the max.r.sd parameter...0.25 seems to be a good value.
calc.copy.number.metrics <- function(samples, params, geno.list, sample.ids=NULL, max.r.sd=NULL) {
    if (!is.null(sample.ids)) {
        samples <- samples[sample.ids,,drop=FALSE]
    }
    if (!is.null(max.r.sd)) {
        params <- .filter.params(params, max.r.sd)
    }
    sapply(rownames(samples), function(sample) {
        g <- geno.list[[sample]]
        metrics <- .calc.copy.number.metrics(params$auto, g)
        if (samples[sample, "sex"] == "m") {
            rbind(metrics, .calc.copy.number.metrics(params$mX, g), .calc.copy.number.metrics(params$Y, g))
        }
        else if (samples[sample, "sex"] == "f") {
            rbind(metrics, .calc.copy.number.metrics(params$fX, g), .calc.copy.number.metrics(params$Y, g))
        }
        else {
            stop(paste("Invalid sex", samples[sample, "sex"], "for sample", sample))
        }
    }, simplify=FALSE)
}

.filter.params <- function(params, max.r.sd) {
    lapply(params, function(x) {
        if (is.null(x)) {
            return(NULL)
        }
        x[x$sAA <= max.r.sd & (is.na(x$sAB) | x$sAB <= max.r.sd) & x$sBB <= max.r.sd,]
    })
}

# Calculate copy number metrics for a single sample and a subset of SNPs.
.calc.copy.number.metrics <- function(params, geno) {
    if (is.null(params)) {
        return(NULL)
    }
    m <- match(rownames(params), rownames(geno))
    params <- params[!is.na(m),]
    geno <- geno[m[!is.na(m)],]

    valid <- geno$x >= 0 & geno$y >= 0
    params <- params[valid,]
    geno <- geno[valid,]

    theta <- (2 / pi) * atan(geno$y / geno$x)
    theta[is.nan(theta)] <- 0
    r <- geno$x + geno$y

    baf <- rep(NaN, length(theta))
    r.exp <- rep(NaN, length(theta))
    s.exp <- rep(NaN, length(theta))

    # BAF is the linear interpolation of the sample theta value between the two nearest
    # canonical clusters. It is clipped to the range [0,1].
    het.missing <- is.na(params$thetaAB)
    if (any(het.missing)) {
        idx1 <- het.missing & theta <= params$thetaAA
        idx2 <- het.missing & theta > params$thetaAA & theta < params$thetaBB
        idx3 <- het.missing & theta >= params$thetaBB

        baf[idx1] <- 0
        baf[idx2] <- ((theta[idx2] - params$thetaAA[idx2]) / (params$thetaBB[idx2] - params$thetaAA[idx2]))
        baf[idx3] <- 1

        r.exp[idx1] <- params$rAA[idx1]
        r.exp[idx2] <- (params$mrA[idx2] * theta[idx2]) + params$brA[idx2]
        r.exp[idx3] <- params$rBB[idx3]

        s.exp[idx1] <- params$sAA[idx1]
        s.exp[idx2] <- (params$msA[idx2] * theta[idx2]) + params$bsA[idx2]
        s.exp[idx3] <- params$sBB[idx3]
    }

    idx1 <- !het.missing & theta <= params$thetaAA
    idx2 <- !het.missing & theta > params$thetaAA & theta < params$thetaAB
    idx3 <- !het.missing & theta == params$thetaAB
    idx4 <- !het.missing & theta > params$thetaAB & theta < params$thetaBB
    idx5 <- !het.missing & theta >= params$thetaBB

    baf[idx1] <- 0
    baf[idx2] <- 0.5 * ((theta[idx2] - params$thetaAA[idx2]) / (params$thetaAB[idx2] - params$thetaAA[idx2]))
    baf[idx3] <- 0.5
    baf[idx4] <- 0.5 + (0.5 * ((theta[idx4] - params$thetaAB[idx4]) / (params$thetaBB[idx4] - params$thetaAB[idx4])))
    baf[idx5] <- 1

    # LRR is the log transformation of the ratio between observed and expected intensity.
    # Expected intensity is calculated by an interpolation of theta along the line between
    # the centroids of the two nearest genotype clusters.
    r.exp[idx1] <- params$rAA[idx1]
    r.exp[idx2] <- (params$mrA[idx2] * theta[idx2]) + params$brA[idx2]
    r.exp[idx3] <- params$rAB[idx3]
    r.exp[idx4] <- (params$mrB[idx4] * theta[idx4]) + params$brB[idx4]
    r.exp[idx5] <- params$rBB[idx5]

    lrr <- log2(r / r.exp)

    # EXPERIMENTAL: investigating the use of Z-scores rather than/in addition to LRR.
    s.exp[idx1] <- params$sAA[idx1]
    s.exp[idx2] <- (params$msA[idx2] * theta[idx2]) + params$bsA[idx2]
    s.exp[idx3] <- params$sAB[idx3]
    s.exp[idx4] <- (params$msB[idx4] * theta[idx4]) + params$bsB[idx4]
    s.exp[idx5] <- params$sBB[idx5]

    z <- (r - r.exp) / s.exp

    valid <- !(is.nan(baf) | is.nan(r.exp) | is.nan(z))
    data.frame(call=geno[valid, "call"], baf=baf[valid], lrr=lrr[valid], z=z[valid], row.names=rownames(geno)[valid], stringsAsFactors=FALSE)
}
