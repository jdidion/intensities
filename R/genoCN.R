# data has five columns: SNP name, chromosome, position, BAF, LRR
cn <- function(name, cn.data, pbf.file, inbred=TRUE, wd=".", cleanup=FALSE, default.state=ifelse(inbred,2,1), 
		default.cn=2, cnv.only=FALSE, params=NULL, ...) {
	# pBs
	pBs <- read.table(pbf.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
	m <- match(cn.data[,1], pBs[,1])
	cn.data <- cn.data[!is.na(m),]
	pBs <- pBs[m[!is.na(m)],4]
	# estimate initial transition probabilities
	data(init.Para.CNA)
	#data(init.Para.CNV)
	if (inbred) {
		init.Para.CNA$trans.begin[1,] <- init.Para.CNA$trans.begin[1,c(2,1,3:9)]
		#init.Para.CNV$trans.m <- init.Para.CNV$trans.m[,c(2,1,3:6)]
		#init.Para.CNV$trans.begin[1,] <- init.Para.CNV$trans.begin[1,c(2,1,3:6)]
	}
	init.Para.CNA$trans.begin <- matrix(init.Para.CNA$trans.begin[1,], length(unique(cn.data$chrm)), 9, byrow=TRUE)
	if (!is.null(params)) {
		for (n in names(params)) {
			init.Para.CNA[[n]] <- params[[n]]
		}
	}
	# set the output dir
	curwd <- getwd()
	tryCatch({
		setwd(wd)
		# run the algorithm
		if (cnv.only) {
			cnv.only <- rep(TRUE, nrow(cn.data))
		}
		else {
			cnv.only <- NULL
		}
		result <- genoCNA(cn.data$name, cn.data$chrm, cn.data$pos, cn.data$lrr, cn.data$baf, pBs, name, init.Para.CNA, outputSeg=TRUE, 
			outputSNP=0, contamination=TRUE, traceIt=0, cnv.only=cnv.only, ...)
		result$name <- name
		#result <- genoCNV(cn.data$name, cn.data$chrm, cn.data$pos, cn.data$lrr, cn.data$baf, pBs, name, init.Para.CNV, outputSeg=TRUE, 
		#	outputSNP=0, loh=TRUE, traceIt=0, ...)
		# read and return the results
		outfile <- paste(name, "segment.txt", sep="_")
		result$segments <- read.table(outfile, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, comment.char='')
		if (is.null(default.cn)) {
			result$ivls <- result$segments[,1:5]
		}
		else {
			result$ivls <- fill.in.missing.intervals(cn.data, result$segments, default.state, default.cn)
		}
		if (cleanup) {
			unlink(outfile)
		}
		result
	}, finally=setwd(curwd))
}

fill.in <- function(sampleIDs, cn.results, default.state, default.cn) {
	names <- as.character(sampleIDs)
	for (n in names) {
		cn.results$cn[[n]]$ivls <- fill.in.missing.intervals(cn.results$data[[n]], cn.results$cn[[n]]$segments, default.state, default.cn)
	}
	cn.results
}

fill.in.missing.intervals <- function(cn.data, segments, default.state, default.cn) {
	ivls <- NULL
	for (chrm in unique(cn.data$chrm)) {
		chrm.ends <- range(cn.data$pos[cn.data$chrm==chrm])
		w <- segments$chr==chrm
		if (any(w)) {
			iv <- segments[w, c(1:5)]
			chrm.ivls <- iv
			n <- nrow(iv)
			if (iv[1,2] > chrm.ends[1]) {
				chrm.ivls <- rbind(chrm.ivls, data.frame(chr=chrm, start=chrm.ends[1], end=iv[1,2]-1, state=default.state, cn=default.cn))
			}
			if (iv[n,3] < chrm.ends[2]) {
				chrm.ivls <- rbind(chrm.ivls, data.frame(chr=chrm, start=iv[n,3]+1, end=chrm.ends[2], state=default.state, cn=default.cn))
			}
			if (n > 1) {
				for (i in 1:(n-1)) {
					if (iv[i,3] < iv[i+1,2]-1) {
						chrm.ivls <- rbind(chrm.ivls, data.frame(chr=chrm, start=iv[i,3]+1, end=iv[i+1,2]-1, state=default.state, cn=default.cn))
					}
				}
			}
		} 
		else {
			chrm.ivls <- data.frame(chr=chrm, start=chrm.ends[1], end=chrm.ends[2], state=default.state, cn=default.cn)
		}
		if (nrow(chrm.ivls) > 1) {
			chrm.ivls <- compress.intervals(chrm.ivls[order(chrm.ivls$start),], val.cols=c(4,5))
		}
		ivls <- rbind(ivls, chrm.ivls)
	}
	ivls
}

summarize <- function(result) {
	summary <- do.call(rbind, by(result$ivls, result$ivls$chr, function(x) tapply(x$end - x$start + 1, x$cn, sum)[as.character(0:4)]))
	summary[is.na(summary)] <- 0
	colnames(summary) <- 0:4
	summary / apply(summary, 1, sum)
}

mean.cn <- function(result) {
	v <- as.vector(by(result$ivls, result$ivls$chr, function(x) {
		sizes <- x$end - x$start + 1
		sum(sizes * x$cn) / sum(sizes)
	}))
	names(v) <- unique(result$ivls$chr)
	v
}

plot.cn.results <- function(sampleIDs, prefix, cn.results, outdir, ...) {
	names <- as.character(sampleIDs)
	for (n in names) {
		jpeg(file.path(outdir, paste(prefix, n, ".jpg", sep="")), width=16, height=12, res=300, units="in")
		plot.cn(cn.results$data[[n]], cn.results$cn[[n]], ...)
		dev.off()
	}
}

plot.cn <- function(cn.data, result, ...) {
	state.colors <- brewer.pal(9, "Set1")
	plot.cn.segments <- function(data, segs, chromosomes, xlim) {
		segs <- segs[segs[,1] %in% names(chromosomes),]
		rng <- range(segs[,5])
		ylim <- c(rng[1] - 0.25, rng[2] + 0.25)
		plot(0:1, type="n", xlim=xlim, ylim=ylim, ylab="Copy Number", xaxt="n", yaxt="n")
		axis(2, seq(rng[1],rng[2]))
		segments(chromosomes[-1], ylim[1], chromosomes[-1], ylim[2], lty=2)
		start <- chromosomes[segs[,1]] + segs[,2]
		end <- chromosomes[segs[,1]] + segs[,3]
		alphas <- as.integer(segs[,1]) %% 2 == 0
		alphas[alphas==0] <- 0.5
		rect(start, segs[,5]-0.25, end, segs[,5]+0.25, col=rgb.alpha(state.colors[segs[,4]], alphas), border=NA)
	}
	plot.sample.copy.number.data(result$name, cn.data, other.plots=list(plot.cn.segments), other.data=result$ivls)
}

write.segments <- function(sampleIDs, results, outfile) {
	write.table(
		do.call(rbind, lapply(as.character(sampleIDs), function(id) cbind(id=rep(id, nrow(results$cn[[id]]$segments)), results$cn[[id]]$segments))),
		outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

est.purity <- function(result) {
	Beta.A.AB = result$mu.b[4,2]
	Beta.B.AB = result$mu.b[4,3]
	 
	# the BAF for genotype AB
	Beta.AB = result$mu.b[1,2]
	 
	# Update BAF estimates to take into account of systematic dye bias
	# this step is not necessary if there is no systematic bias
	Beta.A.AB = 0.5*Beta.A.AB/Beta.AB
	Beta.B.AB = 0.5 + 0.5*(Beta.B.AB - Beta.AB)/(1 - Beta.AB)
	 
	# Re-estimate Beta.A.AB by averaging Beta.A.AB and 1 - Beta.B.AB
	Beta1 = 0.5*(Beta.A.AB + 1 - Beta.B.AB)
	 
	# estimate tumor purity, since we examine the genotype A, nA = 1 and nB = 0
	nA = 1
	nB = 0
	(1 - 2*Beta1)/(Beta1*(nA + nB - 2) + 1 - nB)
}