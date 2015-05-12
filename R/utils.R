# Build a distance matrix from the additive pairwise intensity distances at mulitple SNPs.
# intensities is a four-column data frame of SNP ID, sample ID, x-y, x+y. It is expected to
# be sorted by SNP ID and then sample ID. xform can be used to specify a transformation on
# the distance matricies prior to summation (e.g. log). The result is of class dist.
intensity.dist <- function(intensities, xform=NULL, ...) {
    lapply(sum, by(intensities, intensities[,1], function(x) {
        p <- as.dist(pairdist.default(x[,3] / x[,4], x[,4], ...))
        if (!is.null(xform)) {
            p <- xform(p)
        }
        p
    }))
}

zscore <- function(m) (m - mean(m)) / sd(m)