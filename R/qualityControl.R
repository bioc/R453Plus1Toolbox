##############################
### Read length statistics ###
##############################
.readLengthStats <- function(object) {
  readLengths = width(object)
  stats = vector(mode="numeric", length=5)
  stats[1] = mean(readLengths)
  stats[2] = median(readLengths)
  stats[3] = min(readLengths)
  stats[4] = max(readLengths)
  stats[5] = sd(readLengths)
  stats = round(stats, digits=2)
  names(stats) = c("mean", "median", "min", "max", "sd")
  return(stats)
}

.readLengthStats_sr <- function(object) {
  return(.readLengthStats(sread(object)))
}

.readLengthStats_sff <- function(object) {
  return(.readLengthStats(reads(object)))
}

setMethod("readLengthStats", signature=signature(object="DNAStringSet"), .readLengthStats)
setMethod("readLengthStats", signature=signature(object="ShortRead"), .readLengthStats_sr)
setMethod("readLengthStats", signature=signature(object="SFFContainer"), .readLengthStats_sff)


.readLengthHist <- function(object, cutoff, xlab, ylab, col, breaks, ...) {
  args = list(...)
  readLengths = width(object)
  cut = quantile(readLengths, cutoff)
  readLengths = readLengths[readLengths <= cut]
  par(xaxs = "i")
  par(yaxs = "i")
  if (!is.element("main", names(args))) {
    main = "Histogram of the lengths of the sequences"
    hist(readLengths, main=main, xlab=xlab, ylab=ylab, xlim=c(0, max(readLengths)), breaks=breaks, 
    col=col, ...)
  } else {
    hist(readLengths, xlab=xlab, ylab=ylab, xlim=c(0, max(readLengths)), breaks=breaks, col=col, ...)
  }
}

.readLengthHist_sr <- function(object, cutoff, xlab, ylab, col, breaks, ...) {
  .readLengthHist(sread(object), cutoff, xlab=xlab, ylab=ylab, col=col, breaks=breaks, ...)
}

.readLengthHist_sff <- function(object, cutoff, xlab, ylab, col, breaks, ...) {
  args = list(...)
  if (!is.element("main", names(args))) {
    main = paste("Histogram of the lengths of the sequences \n", name(object), sep="")
    .readLengthHist(reads(object), cutoff, main=main, xlab=xlab, ylab=ylab, col=col, breaks=breaks, ...)
  } else {
    .readLengthHist(reads(object), cutoff, xlab=xlab, ylab=ylab, col=col, breaks=breaks, ...)
  }
}

setMethod("readLengthHist", signature(object="DNAStringSet"), .readLengthHist)
setMethod("readLengthHist", signature(object="ShortRead"), .readLengthHist_sr)
setMethod("readLengthHist", signature(object="SFFContainer"), .readLengthHist_sff)


###############################
### Base quality statistics ###
###############################
.baseQualityStats <- function(object) {
  qualityScores = sapply(as.character(quality(object)), function(x) as.integer(charToRaw(x))-33)
  stats = vector(mode="numeric", length=4)
  stats[1] = sum(sapply(qualityScores, sum)) / sum(width(object))
  stats[2] = min(sapply(qualityScores, min))
  stats[3] = max(sapply(qualityScores, max))
  stats[4] = sqrt(
    sum(sapply(qualityScores, function(x){sum((x-stats[1])^2)})) / 
    (sum(width(object))-1)
    )
  stats = round(stats, digits=2)
  names(stats) = c("mean", "min", "max", "sd")
  return(stats)
}

.baseQualityStats_sr <- function(object) {
  nu = as(quality(object), "numeric")
  stats = vector(mode="numeric", length=4)
  stats[1] = mean(nu)
  stats[2] = min(nu)
  stats[3] = max(nu)
  stats[4] = sd(nu)
  stats = round(stats, digits=2)
  names(stats) = c("mean", "min", "max", "sd")
  return(stats)

}

.baseQualityStats_sff <- function(object) {
  return(.baseQualityStats(reads(object)))
}

setMethod("baseQualityStats", signature(object="QualityScaledDNAStringSet"), .baseQualityStats)
setMethod("baseQualityStats", signature(object="ShortReadQ"), .baseQualityStats_sr)
setMethod("baseQualityStats", signature(object="SFFContainer"), .baseQualityStats_sff)


# Overall quality histogram
.baseQualityHist <- function(object, xlab, ylab, col, breaks, ...) {
  args = list(...)
  onestring = paste(as.character(quality(object)), collapse="")
  qualityScores = as.integer(charToRaw(onestring))-33
  par(xaxs = "i")
  par(yaxs = "i")
  if (!is.element("main", names(args))) {
    main = "Histogram of base quality scores"
    hist(qualityScores, main=main, xlab=xlab, ylab=ylab, breaks=breaks, col=col, ...)
  } else {
    hist(qualityScores, xlab=xlab, ylab=ylab, breaks=breaks, col=col, ...)
  }
}

.baseQualityHist_sr <- function(object, xlab, ylab, col, breaks, ...) {
  .baseQualityHist(as(object, "QualityScaledDNAStringSet"), xlab=xlab, ylab=ylab, col=col, 
    breaks=breaks, ...)
}

.baseQualityHist_sff <- function(object, xlab, ylab, col, breaks, ...) {
  args = list(...)
  if (!is.element("main", names(args))) {
    main = paste("Histogram of base quality scores\n", name(object), sep="")
    .baseQualityHist(reads(object), main=main, xlab=xlab, ylab=ylab, col=col, breaks=breaks, ...)
  } else {
    .baseQualityHist(reads(object), xlab=xlab, ylab=ylab, col=col, breaks=breaks, ...)
  }
}

setMethod("baseQualityHist", signature(object="QualityScaledDNAStringSet"), .baseQualityHist)
setMethod("baseQualityHist", signature(object="ShortReadQ"), .baseQualityHist_sr)
setMethod("baseQualityHist", signature(object="SFFContainer"), .baseQualityHist_sff)


#mean quality per sequence histogram
.sequenceQualityHist <- function(object, xlab, ylab, col, ...) {
  args = list(...)
  qualityScores = sapply(as.character(quality(object)), function(x) as.integer(charToRaw(x))-33)
  if (!is.element("main", names(args)) & !is.element("breaks", names(args))) {
    main = "Histogram of mean quality scores per sequence"
    hist(sapply(qualityScores, mean), main=main, xlab=xlab, ylab=ylab, col=col, breaks=40, ...)
  } else if (!is.element("main", names(args)) & is.element("breaks", names(args))) {
    main = "Histogram of mean quality scores per sequence"
    hist(sapply(qualityScores, mean), main=main, xlab=xlab, ylab=ylab, col=col, ...)
  } else if (is.element("main", names(args)) & !is.element("breaks", names(args))) {
    hist(sapply(qualityScores, mean), xlab=xlab, ylab=ylab, col=col, breaks=40, ...)
  } else {
    hist(sapply(qualityScores, mean), xlab=xlab, ylab=ylab, col=col, ...)
  }
}

.sequenceQualityHist_sr <- function(object, xlab, ylab, col, ...) {
  args = list(...)
  ma = as(quality(object), "matrix")
  if (!is.element("main", names(args)) & !is.element("breaks", names(args))) {
    main = "Histogram of mean quality scores per sequence"
    hist(apply(ma, 1, mean), main=main, xlab=xlab, ylab=ylab, col=col, breaks=max(ma), ...)
  } else if (!is.element("main", names(args)) & is.element("breaks", names(args))) {
    main = "Histogram of mean quality scores per sequence"
    hist(apply(ma, 1, mean), main=main, xlab=xlab, ylab=ylab, col=col, ...)
  } else if (is.element("main", names(args)) & !is.element("breaks", names(args))) {
    hist(apply(ma, 1, mean), xlab=xlab, ylab=ylab, col=col, breaks=max(ma), ...)
  } else {
    hist(apply(ma, 1, mean), xlab=xlab, ylab=ylab, col=col, ...)
  }
}

.sequenceQualityHist_sff <- function(object, xlab, ylab, col, ...) {
  args = list(...)
  if (!is.element("main", names(args))) {
    main = paste("Histogram of mean quality scores per sequence\n", name(object), sep="")
    .sequenceQualityHist(reads(object), main=main, xlab=xlab, ylab=ylab, col=col, ...)
  } else {
    .sequenceQualityHist(reads(object), xlab=xlab, ylab=ylab, col=col, ...)
  } 
}

setMethod("sequenceQualityHist", signature(object="QualityScaledDNAStringSet"), .sequenceQualityHist)
setMethod("sequenceQualityHist", signature(object="ShortReadQ"), .sequenceQualityHist_sr)
setMethod("sequenceQualityHist", signature(object="SFFContainer"), .sequenceQualityHist_sff)


#boxplot of quality per position
.positionQualityBoxplot <- function(object, range, binsize, xlab, ylab, col, ...) {
  args = list(...)
  qs = sapply(as.character(quality(object)), function(x) as.integer(charToRaw(x))-33)
  nreads = length(object)
  readlengths = width(object)
  nbases = max(readlengths)
  ma = matrix(NA, nrow=nreads, ncol=nbases) # empty matrix, size based on longest read
  for(i in 1:nreads) {
    ma[i, seq(along=qs[[i]])] = qs[[i]]
  }
  if (missing(range)) {
    range = c(1, nbases)
  }
  if(length(range) == 1) {
    start = 1
    end = range
  } else {
    start = range[1]
    end = range[2]
  }
  if(start > dim(ma)[2] | end > dim(ma)[2]) {
    stop("Values are out of range")
  }
  ma = ma[, start:end]
  goalsize = binsize * ceiling(dim(ma)[2] / binsize) # size dividable by binsize
  mb = matrix(NA, nrow=nreads, ncol=goalsize-dim(ma)[2]) # fill with NAs
  mc = cbind(ma, mb)
  md = matrix(mc, ncol=goalsize/binsize)
  if (!is.element("main", names(args))) {
    main = "Boxplot of quality per position"
    boxplot(md, outline=FALSE, col=col, axes=FALSE, main=main, xlab=xlab, ylab=ylab, frame=TRUE, 
    ylim=c(0, max(md, na.rm=TRUE)), ...)
  } else {
    boxplot(md, outline=FALSE, col=col, axes=FALSE, xlab=xlab, ylab=ylab, frame=TRUE, 
    ylim=c(0, max(md, na.rm=TRUE)), ...)
  }
  myLabels <- pretty(0:dim(ma)[2], n=4, min.n=1)
  myTicks <- myLabels/binsize
  axis(1, at=myTicks, labels=myLabels)
  axis(2, at=c(0, 10, 20, 30, 40), labels=TRUE)
}

.positionQualityBoxplot_sr <- function(object, range, binsize, xlab, ylab, col, ...) {
  args = list(...)
  ma = as(quality(object), "matrix")
  if (missing(range)) {
    range = c(1, dim(ma)[2])
  }
  if(length(range) == 1) {
    start = 1
    end = range
  } else {
    start = range[1]
    end = range[2]
  }
  if(start > dim(ma)[2] | end > dim(ma)[2]) {
    stop("Values are out of range")
  }
  ma = ma[, start:end]
  goalsize = binsize * ceiling(dim(ma)[2] / binsize) # size dividable by binsize
  mb = matrix(NA, nrow=dim(ma)[1], ncol=goalsize-dim(ma)[2]) # fill with NAs
  mc = cbind(ma, mb)
  md = matrix(mc, ncol=goalsize/binsize)
  maxq = max(md, na.rm=TRUE)
  if (!is.element("main", names(args))) {
    main = "Boxplot of quality per position"
    boxplot(md, outline=FALSE, col=col, axes=FALSE, main=main, xlab=xlab, ylab=ylab, frame=TRUE, 
    ylim=c(0, maxq), ...)
  } else {
    boxplot(md, outline=FALSE, col=col, axes=FALSE, xlab=xlab, ylab=ylab, frame=TRUE, 
    ylim=c(0, maxq), ...)
  }  
  myLabels <- pretty(0:dim(ma)[2], n=4, min.n=1)
  myTicks <- myLabels/binsize
  axis(1, at=myTicks, labels=myLabels)
  axis(2, at=c(0, 10, 20, 30, 40), labels=TRUE)
}

.positionQualityBoxplot_sff <- function(object, range, binsize, xlab, ylab, col, ...) {
  args = list(...)
  if (!is.element("main", names(args))) {
    main = paste("Boxplot of quality per position\n", name(object), sep="")
    .positionQualityBoxplot(reads(object), range, binsize, main=main, xlab=xlab, ylab=ylab, col=col, ...)
  } else {
    .positionQualityBoxplot(reads(object), range, binsize, xlab=xlab, ylab=ylab, col=col, ...)
  }
}

setMethod("positionQualityBoxplot", signature(object="QualityScaledDNAStringSet"), .positionQualityBoxplot)
setMethod("positionQualityBoxplot", signature(object="ShortReadQ"), .positionQualityBoxplot_sr)
setMethod("positionQualityBoxplot", signature(object="SFFContainer"), .positionQualityBoxplot_sff)


####################################
### Base distribution statistics ###
####################################
.baseFrequency <- function(object) {
  af1 = alphabetFrequency(object, baseOnly=TRUE, collapse=TRUE)
  af2 = alphabetFrequency(object, as.prob=TRUE, baseOnly=TRUE, collapse=TRUE) * 100
  af1["total"] = sum(af1)
  af2["total"] = 100
  af2 = round(af2, digits=2)
  ff = data.frame(absolute=af1, relative=af2)
  return(ff)
}

.baseFrequency_sr <- function(object) {
  return(.baseFrequency(sread(object)))
}

.baseFrequency_sff <- function(object) {
  re = reads(object)
  return(.baseFrequency(re))
}

setMethod("baseFrequency", signature=signature(object="DNAStringSet"), .baseFrequency)
setMethod("baseFrequency", signature=signature(object="ShortRead"), .baseFrequency_sr)
setMethod("baseFrequency", signature=signature(object="SFFContainer"), .baseFrequency_sff)


.nucleotideCharts <- function(object, range, linetypes, linecols, xlab, ylab, ...) {
  args = list(...)
  readLengths = width(object)
  if(length(range) == 1) {
    start = 1
    end = quantile(readLengths, range)
  } else {
    start = range[1]
    end = range[2]
  }
  abc = alphabetByCycle(object, c("A", "C", "G", "T", "N"))
  nucsPerPos = colSums(abc)
  if(start > dim(abc)[2] | end > dim(abc)[2]) {
    stop("Values are out of range")
  }
  if (!is.element("main", names(args))) {
    main = "Nucleotide frequency"
    plot(x=start:end, y=abc["A", start:end]/nucsPerPos[start:end], type=linetypes["A"], 
      col=linecols["A"], main=main, xlab=xlab, ylab=ylab, ylim=c(0,1), ...)
  } else {
    plot(x=start:end, y=abc["A", start:end]/nucsPerPos[start:end], type=linetypes["A"], 
      col=linecols["A"], xlab=xlab, ylab=ylab, ylim=c(0,1), ...)
  }
  lines(x=start:end, y=abc["C", start:end]/nucsPerPos[start:end], type=linetypes["C"], 
    col=linecols["C"])
  lines(x=start:end, y=abc["G", start:end]/nucsPerPos[start:end], type=linetypes["G"], 
    col=linecols["G"])
  lines(x=start:end, y=abc["T", start:end]/nucsPerPos[start:end], type=linetypes["T"], 
    col=linecols["T"])
  lines(x=start:end, y=abc["N", start:end]/nucsPerPos[start:end], type=linetypes["N"], 
    col=linecols["N"])
  legend(x="top", horiz=TRUE, legend=c("A", "C", "G", "T", "N"), 
    fill=c(linecols["A"], linecols["C"], linecols["G"], linecols["T"], linecols["N"]))
}

.nucleotideCharts_sr <- function(object, range, linetypes, linecols, xlab, ylab, ...) {
  .nucleotideCharts(sread(object), range, linetypes=linetypes, linecols=linecols, xlab=xlab, 
    ylab=ylab, ...)
}

.nucleotideCharts_sff <- function(object, range, linetypes, linecols, xlab, ylab, ...) {
  args = list(...)
  if (!is.element("main", names(args))) {
    main = paste("Nucleotide frequency\n", name(object), sep="")
    .nucleotideCharts(reads(object), range, linetypes=linetypes, linecols=linecols, main=main, 
      xlab=xlab, ylab=ylab, ...)
  } else {
    .nucleotideCharts(reads(object), range, linetypes=linetypes, linecols=linecols,
      xlab=xlab, ylab=ylab, ...)
  } 
}

setMethod("nucleotideCharts", signature(object="DNAStringSet"), .nucleotideCharts)
setMethod("nucleotideCharts", signature(object="ShortRead"), .nucleotideCharts_sr)
setMethod("nucleotideCharts", signature(object="SFFContainer"), .nucleotideCharts_sff)


.gcContent <- function(object) {
  af = alphabetFrequency(object, baseOnly=TRUE, collapse=TRUE)
  gc = ((af["G"] + af["C"]) / (af["A"] + af["C"] + af["G"] + af["T"])) * 100
  names(gc) = "GCcontent"
  gc = round(gc, digits=2)
  return(gc)
}

.gcContent_sr <- function(object) {
  return(.gcContent(sread(object)))
}

.gcContent_sff <- function(object) {
  reads = reads(object)
  return(.gcContent(reads))
}

setMethod("gcContent", signature(object="DNAStringSet"), .gcContent)
setMethod("gcContent", signature(object="ShortRead"), .gcContent_sr)
setMethod("gcContent", signature(object="SFFContainer"), .gcContent_sff)


.gcPerPosition <- function(object, range, type, col, xlab, ylab, ...) {
  args = list(...)
  readLengths = width(object)
  if(length(range) == 1) {
    start = 1
    end = quantile(readLengths, range)
  } else {
    start = range[1]
    end = range[2]
  }
  abc = alphabetByCycle(object, c("A", "C", "G", "T", "N"))
  nucsPerPos = colSums(abc)
  if(start > dim(abc)[2] | end > dim(abc)[2]) {
    stop("Values are out of range")
  }
  if (!is.element("main", names(args))) {
    main = "GC content per base"
    plot(x=start:end, y=(abc["G", start:end] + abc["C", start:end])/nucsPerPos[start:end], 
      type=type, col=col, main=main, xlab=xlab, ylab=ylab, ylim=c(0,1), ...)
  } else {
    plot(x=start:end, y=(abc["G", start:end] + abc["C", start:end])/nucsPerPos[start:end], 
      type=type, col=col, xlab=xlab, ylab=ylab, ylim=c(0,1), ...)
  }
}

.gcPerPosition_sr <- function(object, range, type, col, xlab, ylab, ...) {
  .gcPerPosition(sread(object), range, type=type, col=col, xlab=xlab, ylab=ylab, ...)
}

.gcPerPosition_sff <- function(object, range, type, col, xlab, ylab, ...) {
  args = list(...)
  if (!is.element("main", names(args))) {
    main = paste("GC content per base\n", name(object), sep="")
    .gcPerPosition(reads(object), range, main=main, type=type, col=col, xlab=xlab, ylab=ylab, ...)
  } else {
    .gcPerPosition(reads(object), range, type=type, col=col, xlab=xlab, ylab=ylab, ...)
  }
}

setMethod("gcPerPosition", signature(object="DNAStringSet"), .gcPerPosition)
setMethod("gcPerPosition", signature(object="ShortRead"), .gcPerPosition_sr)
setMethod("gcPerPosition", signature(object="SFFContainer"), .gcPerPosition_sff)


.gcContentHist <- function(object, xlab, ylab, col, breaks, ...) {
  args = list(...)
  af = alphabetFrequency(object, baseOnly=TRUE)
  gc = ((af[,"G"] + af[,"C"]) / (af[,"A"] + af[,"C"] + af[,"G"] + af[,"T"])) * 100
  if (!is.element("main", names(args))) {
    main = "Histogram of GC content per sequence"
    hist(gc, main=main, xlim=c(0,100), xlab=xlab, ylab=ylab, col=col, breaks=breaks, ...)
  } else {
    hist(gc, xlim=c(0,100), xlab=xlab, ylab=ylab, col=col, breaks=breaks, ...)
  }
}

.gcContentHist_sr <- function(object, xlab, ylab, col, breaks, ...) {
  .gcContentHist(sread(object), xlab=xlab, ylab=ylab, col=col, breaks=breaks, ...)
}

.gcContentHist_sff <- function(object, xlab, ylab, col, breaks, ...) {
  args = list(...)
  if (!is.element("main", names(args))) {
    main = paste("Histogram of GC content per sequence\n", name(object), sep="")
    .gcContentHist(reads(object), main=main, xlab=xlab, ylab=ylab, col=col, breaks=breaks, ...)
  } else {
    .gcContentHist(reads(object), xlab=xlab, ylab=ylab, col=col, breaks=breaks, ...)
  }
}

setMethod("gcContentHist", signature(object="DNAStringSet"), .gcContentHist)
setMethod("gcContentHist", signature(object="ShortRead"), .gcContentHist_sr)
setMethod("gcContentHist", signature(object="SFFContainer"), .gcContentHist_sff)


###########################
### Sequence complexity ###
###########################

### DUST approach ###
.complexity.dust <- function(object, xlab, ylab, xlim, col, breaks, ...) {
  args = list(...)
  scores = vector(mode="numeric", length=length(object))
  tfq = trinucleotideFrequency(object) # matrix mit nreads rows and 64 cols
  rls = width(object)
  scaling = ((rls - 2) * (rls - 2 - 1)) / ((rls - 2 - 1) * 2)
  scores = rowSums(tfq * (tfq - 1)) * 100 / (2 * scaling * (rls - 2 - 1))
  if (!is.element("main", names(args))) {
    main = "Histogram of DUST complexity scores"
    hist(scores, main=main, xlab=xlab, ylab=ylab, xlim=xlim, col=col, breaks=breaks, ...)
  } else {
    hist(scores, xlab=xlab, ylab=ylab, xlim=xlim, col=col, breaks=breaks, ...)
  }
  abline(v=7, col="blue")
  return(scores)
}

.complexity.dust_sr <- function(object, xlab, ylab, xlim, col, breaks, ...) {
  .complexity.dust(sread(object), xlab=xlab, ylab=ylab, xlim=xlim, col=col, breaks=breaks, ...)
}

.complexity.dust_sff <- function(object, xlab, ylab, xlim, col, breaks, ...) {
  args = list(...)
  if (!is.element("main", names(args))) {
    main = paste("Histogram of DUST complexity scores\n", name(object), sep="")
    .complexity.dust(reads(object), main=main, xlab=xlab, ylab=ylab, xlim=xlim, col=col, 
      breaks=breaks, ...)
  } else {
    .complexity.dust(reads(object), xlab=xlab, ylab=ylab, xlim=xlim, col=col, breaks=breaks, ...)
  }
}

setMethod("complexity.dust", signature(object="DNAStringSet"), .complexity.dust)
setMethod("complexity.dust", signature(object="ShortRead"), .complexity.dust_sr)
setMethod("complexity.dust", signature(object="SFFContainer"), .complexity.dust_sff)

### Entropy approach ###
.complexity.entropy <- function(object, xlab, ylab, xlim, col, breaks, ...) {
  args = list(...)
  scores = vector(mode="numeric", length=length(object))
  tfq = trinucleotideFrequency(object) # matrix with nreads rows and 64 cols
  rls = width(object)
  fac1 = tfq / (rls - 2)
  fac2 = ifelse(fac1 == 0, 0, log(fac1, base=ifelse(rls < 66, rls - 2, 64)))
  scores = (-100) * rowSums(fac1 * fac2)
  if (!is.element("main", names(args))) {
    main = "Histogram of Entropy complexity scores"
    hist(scores, main=main, xlab=xlab, ylab=ylab, xlim=xlim, breaks=breaks, col=col, ...)
  } else {
    hist(scores, xlab=xlab, ylab=ylab, xlim=xlim, breaks=breaks, col=col, ...)
  }
  abline(v=70, col="blue")
  return(scores)
}

.complexity.entropy_sr <- function(object, xlab, ylab, xlim, col, breaks, ...) {
  .complexity.entropy(sread(object), xlab=xlab, ylab=ylab, xlim=xlim, col=col, breaks=breaks, ...)
}

.complexity.entropy_sff <- function(object, xlab, ylab, xlim, col, breaks, ...) {
  args = list(...)
  if (!is.element("main", names(args))) {
    main = paste("Histogram of Entropy complexity scores\n", name(object), sep="")
    .complexity.entropy(reads(object), main=main, xlab=xlab, ylab=ylab, xlim=xlim, col=col, 
      breaks=breaks, ...)
  } else {
    .complexity.entropy(reads(object), xlab=xlab, ylab=ylab, xlim=xlim, col=col, breaks=breaks, ...)
  }
}

setMethod("complexity.entropy", signature(object="DNAStringSet"), .complexity.entropy)
setMethod("complexity.entropy", signature(object="ShortRead"), .complexity.entropy_sr)
setMethod("complexity.entropy", signature(object="SFFContainer"), .complexity.entropy_sff)

###############################
### Dinucleotide odds ratio ###
###############################
.dinucleotideOddsRatio <- function(object, xlab, col, ...) {
  args = list(...)
  dfq = dinucleotideFrequency(object) # matrix with nreads rows and 16 cols
  nfq = oligonucleotideFrequency(object, width=1) # matrix with nreads rows and 4 cols
  dinucs = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", 
    "TG", "TT")
  rcdinucs = c("TT", "GT", "GC", "AT", "TG", "GG", "CG", "AG", "TC", "GC", "CC", "AC", "TA", "GA", 
    "CA", "AA")
  scorema = matrix(nrow=length(object), ncol=length(dinucs))
  colnames(scorema) = dinucs
  for(i in 1:length(dinucs)) {
    div1 = dfq[, dinucs[i]] + dfq[, rcdinucs[i]]
    div2 = nfq[, substr(dinucs[i], 1, 1)] + nfq[, substr(rcdinucs[i], 2, 2)]
    div3 = nfq[, substr(dinucs[i], 2, 2)] + nfq[, substr(rcdinucs[i], 1, 1)]
    div4 = 2 * (rowSums(nfq)^2) / rowSums(dfq)
    scorema[, dinucs[i]] = (div1 / (div2 * div3)) * div4
  }
  means = apply(scorema, 2, mean)
  if (!is.element("main", names(args))) {
    main = "Dinucleotide odds ratio"
    bar = barplot(means-1, main=main, xlab=xlab, horiz=TRUE, xlim=c(-1,1), col=col, las=1, 
      axes=FALSE, axisnames=FALSE, ...)
  } else {
    bar = barplot(means-1, xlab=xlab, horiz=TRUE, xlim=c(-1,1), col=col, las=1, axes=FALSE, 
      axisnames=FALSE, ...)
  }
  abline(v=0, col="black")
  axis(1, at=pretty(-1:1, n=4, min.n=1), labels=pretty(0:2, n=4, min.n=1))
  text(x=means-1+ifelse(means < 1, -0.08, 0.08), y=bar[, 1], labels=dinucs)
  return(scorema)
}

.dinucleotideOddsRatio_sr <- function(object, xlab, col, ...) {
  .dinucleotideOddsRatio(sread(object), xlab=xlab, col=col, ...)
}

.dinucleotideOddsRatio_sff <- function(object, xlab, col, ...) {
  args = list(...)
  if (!is.element("main", names(args))) {
    main = paste("Dinucleotide odds ratio\n", name(object), sep="")
    .dinucleotideOddsRatio(reads(object), main=main, xlab=xlab, col=col, ...)
  } else {
    .dinucleotideOddsRatio(reads(object), xlab=xlab, col=col, ...)
  }
}

setMethod("dinucleotideOddsRatio", signature(object="DNAStringSet"), .dinucleotideOddsRatio)
setMethod("dinucleotideOddsRatio", signature(object="ShortRead"), .dinucleotideOddsRatio_sr)
setMethod("dinucleotideOddsRatio", signature(object="SFFContainer"), .dinucleotideOddsRatio_sff)


###########################
### Flowgram statistics ###
###########################

.flowgramBarplot <- function(x, range, xlab, ylab, col, ...) {
  args = list(...)
  nflows = length(flowgram(x))
  if(length(range) != 2 | range[1] < 0 | range[2] > nflows | range[1] > range[2]) {
    stop("Illegal range values.")
  }
  flowChars = unlist(strsplit(flowChars(x), NULL))[range[1]:range[2]]
  flowgram = flowgram(x)[range[1]:range[2]]
  colors = col[match(flowChars, names(col))]
  if (!is.element("main", names(args))) {
    main = paste("Flowgram intensity\n", name(x), sep="")
    barplot(flowgram, main=main, xlab=xlab, ylab=ylab, col=colors, border=colors, width=1, space=0, ...)
  } else {
    barplot(flowgram, xlab=xlab, ylab=ylab, col=colors, border=colors, width=1, space=0, ...)
  }
  myLabels <- pretty(range[1]:range[2], n=5)
  myTicks <- myLabels-range[1]
  axis(1, at=myTicks, labels=myLabels)
  legend(x="top", horiz=TRUE, legend=c("A", "C", "G", "T"), 
    fill=c(col["A"], col["C"], col["G"], col["T"]))
}

setMethod("flowgramBarplot", signature(x="SFFRead"), .flowgramBarplot)


.homopolymerHist <- function(x, range, xlab, ylab, col, ...) {
  args = list(...)
  flowChars = unlist(strsplit(flowChars(x), NULL))[range[1]:range[2]]
  flowgram = round(flowgram(x)[range[1]:range[2]]/100)
  longest = max(flowgram)
  ta = table(factor(flowgram[flowChars == "A"], levels=0:longest))
  tc = table(factor(flowgram[flowChars == "C"], levels=0:longest))
  tg = table(factor(flowgram[flowChars == "G"], levels=0:longest))
  tt = table(factor(flowgram[flowChars == "T"], levels=0:longest))
  ma = matrix(c(ta,tc,tg,tt), nrow=4, byrow=TRUE)
  rownames(ma) = c("A", "C", "G", "T")
  colnames(ma) = 0:longest
  if (!is.element("main", names(args))) {
    main = paste("Homopolymers in read:\n", name(x), sep="")
    barplot(ma, beside=TRUE, main=main, xlab=xlab, ylab=ylab, col=col, ...)
  } else {
    barplot(ma, beside=TRUE, xlab=xlab, ylab=ylab, col=col, ...)
  }
  legend(x="topright", legend=c("A", "C", "G", "T"), 
    fill=c(col["A"], col["C"], col["G"], col["T"]))
  #as.data.frame(table(flowgram[flowChars == "A"], exclude=c(0,1)))
}

setMethod("homopolymerHist", signature(x="SFFRead"), .homopolymerHist)
