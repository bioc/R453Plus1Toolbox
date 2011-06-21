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


.readLengthHist <- function(object, cutoff, name) {
  if (missing(name)) {
    name = ""
  }
  readLengths = width(object)
  cut = quantile(readLengths, cutoff)
  readLengths = readLengths[readLengths <= cut]
  par(xaxs = "i")
  par(yaxs = "i")
  title = paste("Histogram of the lengths of the sequences \n", name, sep="")
  hist(readLengths, main=title, xlab="Read length", 
    ylab="Number of sequences", xlim=c(0, max(readLengths)), breaks=100, col="firebrick1")
}

.readLengthHist_sr <- function(object, cutoff, name) {
  .readLengthHist(sread(object), cutoff, name)
}

.readLengthHist_sff <- function(object, cutoff, name) {
  if (missing(name)) {
    name = filename(object)
  }
  reads = reads(object)
  .readLengthHist(reads, cutoff, name)
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
.baseQualityHist <- function(object, name) {
  if (missing(name)) {
    name = ""
  }
  onestring = paste(as.character(quality(object)), collapse="")
  qualityScores = as.integer(charToRaw(onestring))-33
  par(xaxs = "i")
  par(yaxs = "i")
  title = paste("Histogram of base quality scores\n", name, sep="")
  hist(qualityScores, main=title, xlab="Quality score", 
    ylab="Number of bases", breaks=40, col="firebrick1")
}

.baseQualityHist_sr <- function(object, name) {
  .baseQualityHist(as(object, "QualityScaledDNAStringSet"), name)
}

.baseQualityHist_sff <- function(object, name) {
  if (missing(name)) {
    name = filename(object)
  }
  reads = reads(object)
  .baseQualityHist(reads, name)
}

setMethod("baseQualityHist", signature(object="QualityScaledDNAStringSet"), .baseQualityHist)
setMethod("baseQualityHist", signature(object="ShortReadQ"), .baseQualityHist_sr)
setMethod("baseQualityHist", signature(object="SFFContainer"), .baseQualityHist_sff)


#mean quality per sequence histogram
.sequenceQualityHist <- function(object, name) {
  if (missing(name)) {
    name = ""
  }
  qualityScores = sapply(as.character(quality(object)), function(x) as.integer(charToRaw(x))-33)
  title = paste("Histogram of mean quality scores per sequence\n", name, sep="")
  hist(sapply(qualityScores, mean), 
    main=title, 
    xlab="Mean of quality scores per sequence", ylab="Number of sequences", col="firebrick1", 
    breaks=40)
}

.sequenceQualityHist_sr <- function(object, name) {
  if (missing(name)) {
    name = ""
  }
  ma = as(quality(object), "matrix")
  title = paste("Histogram of mean quality scores per sequence\n", name, sep="")
  hist(apply(ma, 1, mean), 
    main=title, 
    xlab="Mean of quality scores per sequence", ylab="Number of sequences", col="firebrick1", 
    breaks=max(ma))
}

.sequenceQualityHist_sff <- function(object, name) {
  if (missing(name)) {
    name = filename(object)
  }
  reads = reads(object)
  .sequenceQualityHist(reads, name)
}

setMethod("sequenceQualityHist", signature(object="QualityScaledDNAStringSet"), .sequenceQualityHist)
setMethod("sequenceQualityHist", signature(object="ShortReadQ"), .sequenceQualityHist_sr)
setMethod("sequenceQualityHist", signature(object="SFFContainer"), .sequenceQualityHist_sff)


#boxplot of quality per position
.positionQualityBoxplot <- function(object, range, binsize, name) {
  if (missing(name)) {
    name = ""
  }
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
  title = paste("Boxplot of quality per position\n", name, sep="")
  boxplot(md, outline=FALSE, col="firebrick1", axes=FALSE, main=title, 
    xlab=paste("Read position in bp (Bin size: ", binsize, "bp)", sep=""), ylab="Quality score", 
    frame=TRUE, ylim=c(0, max(md, na.rm=TRUE)))
  myLabels <- pretty(0:dim(ma)[2], n=4, min.n=1)
  myTicks <- myLabels/binsize
  axis(1, at=myTicks, labels=myLabels)
  axis(2, at=c(0, 10, 20, 30, 40), labels=TRUE)
}

.positionQualityBoxplot_sr <- function(object, range, binsize, name) {
  if (missing(name)) {
    name = ""
  }
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
  title = paste("Boxplot of quality per position\n", name, sep="")
  boxplot(md, outline=FALSE, col="firebrick1", axes=FALSE, main=title, 
    xlab=paste("Read position in bp (Bin size: ", binsize, "bp)", sep=""), ylab="Quality score", 
    frame=TRUE, ylim=c(0, maxq))    
  myLabels <- pretty(0:dim(ma)[2], n=4, min.n=1)
  myTicks <- myLabels/binsize
  axis(1, at=myTicks, labels=myLabels)
  axis(2, at=c(0, 10, 20, 30, 40), labels=TRUE)
}

.positionQualityBoxplot_sff <- function(object, range, binsize, name) {
  if (missing(name)) {
    name = filename(object)
  }
  reads = reads(object)
  .positionQualityBoxplot(reads, range, binsize, name)
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


.nucleotideCharts <- function(object, range, name) {
  if (missing(name)) {
    name = ""
  }
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
  title = paste("Nucleotide frequency\n", name, sep="")
  plot(x=start:end, y=abc["A", start:end]/nucsPerPos[start:end], type="l", col="black",
    main=title, xlab="Position", ylab="Frequency", ylim=c(0,1))
  lines(x=start:end, y=abc["C", start:end]/nucsPerPos[start:end], type="l", col="red")
  lines(x=start:end, y=abc["G", start:end]/nucsPerPos[start:end], type="l", col="blue")
  lines(x=start:end, y=abc["T", start:end]/nucsPerPos[start:end], type="l", col="green")
  lines(x=start:end, y=abc["N", start:end]/nucsPerPos[start:end], type="l", col="grey70")
  legend(x="top", horiz=TRUE, legend=c("A", "C", "G", "T", "N"), fill=c("black", "red", "blue", "green", 
    "grey70"))
}

.nucleotideCharts_sr <- function(object, range, name) {
  .nucleotideCharts(sread(object), range, name)
}

.nucleotideCharts_sff <- function(object, range, name) {
  if (missing(name)) {
    name = filename(object)
  }
  reads = reads(object)
  .nucleotideCharts(reads, range, name)
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


.gcPerPosition <- function(object, range, name) {
  if (missing(name)) {
    name = ""
  }
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
  title = paste("GC content per base\n", name, sep="")
  plot(x=start:end, y=(abc["G", start:end] + abc["C", start:end])/nucsPerPos[start:end], type="l", 
    col="blue", main=title, xlab="Position", ylab="Frequency", ylim=c(0,1))
}

.gcPerPosition_sr <- function(object, range, name) {
  .gcPerPosition(sread(object), range, name)
}

.gcPerPosition_sff <- function(object, range, name) {
  if (missing(name)) {
    name = filename(object)
  }
  reads = reads(object)
  .gcPerPosition(reads, range, name)
}

setMethod("gcPerPosition", signature(object="DNAStringSet"), .gcPerPosition)
setMethod("gcPerPosition", signature(object="ShortRead"), .gcPerPosition_sr)
setMethod("gcPerPosition", signature(object="SFFContainer"), .gcPerPosition_sff)


.gcContentHist <- function(object, name) {
  if (missing(name)) {
    name = ""
  }
  af = alphabetFrequency(object, baseOnly=TRUE)
  gc = ((af[,"G"] + af[,"C"]) / (af[,"A"] + af[,"C"] + af[,"G"] + af[,"T"])) * 100
  title = paste("Histogram of GC content per sequence\n", name, sep="")
  hist(gc, main=title, xlim=c(0,100),
    xlab="GC content", ylab="Number of sequences", col="firebrick1", 
    breaks=50)
}

.gcContentHist_sr <- function(object, name) {
  .gcContentHist(sread(object), name)
}

.gcContentHist_sff <- function(object, name) {
  if (missing(name)) {
    name = filename(object)
  }
  reads = reads(object)
  .gcContentHist(reads, name)
}

setMethod("gcContentHist", signature(object="DNAStringSet"), .gcContentHist)
setMethod("gcContentHist", signature(object="ShortRead"), .gcContentHist_sr)
setMethod("gcContentHist", signature(object="SFFContainer"), .gcContentHist_sff)


###########################
### Sequence complexity ###
###########################

### DUST approach ###
.complexity.dust <- function(object, name) {
  if (missing(name)) {
    name = ""
  }
  scores = vector(mode="numeric", length=length(object))
  tfq = trinucleotideFrequency(object) # matrix mit nreads rows and 64 cols
  rls = width(object)
  scaling = ((rls - 2) * (rls - 2 - 1)) / ((rls - 2 - 1) * 2)
  scores = rowSums(tfq * (tfq - 1)) * 100 / (2 * scaling * (rls - 2 - 1))
  title = paste("Histogram of DUST complexity scores\n", name, sep="")
  hist(scores, main=title, 
    xlab="Complexity score (0=high, 100=low)", ylab="Number of sequences", xlim=c(0, 100), 
    breaks=100, col="firebrick1")
  abline(v=7, col="blue")
  return(scores)
}

.complexity.dust_sr <- function(object, name) {
  .complexity.dust(sread(object), name)
}

.complexity.dust_sff <- function(object, name) {
  if (missing(name)) {
    name = filename(object)
  }
  reads = reads(object)
  .complexity.dust(reads, name)
}

setMethod("complexity.dust", signature(object="DNAStringSet"), .complexity.dust)
setMethod("complexity.dust", signature(object="ShortRead"), .complexity.dust_sr)
setMethod("complexity.dust", signature(object="SFFContainer"), .complexity.dust_sff)

### Entropy approach ###
.complexity.entropy <- function(object, name) {
  if (missing(name)) {
    name = ""
  }
  scores = vector(mode="numeric", length=length(object))
  tfq = trinucleotideFrequency(object) # matrix with nreads rows and 64 cols
  rls = width(object)
  fac1 = tfq / (rls - 2)
  fac2 = ifelse(fac1 == 0, 0, log(fac1, base=ifelse(rls < 66, rls - 2, 64)))
  scores = (-100) * rowSums(fac1 * fac2)
  title = paste("Histogram of Entropy complexity scores\n", name, sep="")
  hist(scores, main=title, 
    xlab="Complexity score (0=low, 100=high)", 
    ylab="Number of sequences", xlim=c(0, 100), breaks=100, col="firebrick1")
  abline(v=70, col="blue")
  return(scores)
}

.complexity.entropy_sr <- function(object, name) {
  .complexity.entropy(sread(object), name)
}

.complexity.entropy_sff <- function(object, name) {
  if (missing(name)) {
    name = filename(object)
  }
  reads = reads(object)
  .complexity.entropy(reads, name)
}

setMethod("complexity.entropy", signature(object="DNAStringSet"), .complexity.entropy)
setMethod("complexity.entropy", signature(object="ShortRead"), .complexity.entropy_sr)
setMethod("complexity.entropy", signature(object="SFFContainer"), .complexity.entropy_sff)

###############################
### Dinucleotide odds ratio ###
###############################
.dinucleotideOddsRatio <- function(object, name) {
  if (missing(name)) {
    name = ""
  }
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
  title = paste("Dinucleotide odds ratio\n", name, sep="")
  bar = barplot(means-1, main=title, 
    xlab="Under-/over-representation of dinucleotides", horiz=TRUE, xlim=c(-1,1), col="firebrick1", 
    las=1, axes=FALSE, axisnames=FALSE)
  abline(v=0, col="black")
  axis(1, at=pretty(-1:1, n=4, min.n=1), labels=pretty(0:2, n=4, min.n=1))
  text(x=means-1+ifelse(means < 1, -0.08, 0.08), y=bar[, 1], labels=dinucs)
  return(scorema)
}

.dinucleotideOddsRatio_sr <- function(object, name) {
  .dinucleotideOddsRatio(sread(object), name)
}

.dinucleotideOddsRatio_sff <- function(object, name) {
  if (missing(name)) {
    name = filename(object)
  }
  reads = reads(object)
  .dinucleotideOddsRatio(reads, name)
}

setMethod("dinucleotideOddsRatio", signature(object="DNAStringSet"), .dinucleotideOddsRatio)
setMethod("dinucleotideOddsRatio", signature(object="ShortRead"), .dinucleotideOddsRatio_sr)
setMethod("dinucleotideOddsRatio", signature(object="SFFContainer"), .dinucleotideOddsRatio_sff)
