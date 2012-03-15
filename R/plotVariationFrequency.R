.plotVariationFrequency_data.frame <- function(object, plotRange, col, sequenceCex, xlab, ylab, ylim, ...) {

    if (missing(sequenceCex)) {
        sequenceCex = 0.5
    }
    if (missing(col)) {
        col = c("green", "blue", "black", "red", "purple", "gray", "blue")
    }
    if (missing(ylab)) {
        ylab = "Variation (%)"
    }
    # xlab default value is defined below, because plotRange has be checked first.
  
    refSeq = object[,2]                        # Refposn.
    refPos = floor(object[,1])                 # X
    subMatrix = as.matrix(object[, 3:8])       # A.. C.., G.., T.., N.., X...
    numReads = object[, 9]                     # X.1
    rm(object)

    # match given ploting range (bp pos in amplicon ref) to row in data.frame
    if (min(refPos) > plotRange[1] | max(refPos) < plotRange[2]) {
       stop("plotRange is out of bounce.")
    }
    plotRange = match(plotRange[1], refPos):match(plotRange[2], refPos)

    if (missing(xlab)) {
        xlab = paste("Reference sequence (position ", refPos[plotRange[1]], "-",
            refPos[plotRange[length(plotRange)]], ")", sep="")
    }


    # define the y-axis in 5% steps
    if (missing(ylim)) {
        maxPercent = ceiling(max(t(subMatrix)[, plotRange])/10 * 2) * 10 / 2
        if (maxPercent > 80) {
            maxPercent = 100
        }
        ylim = c(0, maxPercent)
    } else {
        maxPercent = ylim[2]
    }
 

    par(mar=c(5.1, 4.1, 4.1, 5.1))
    barplot(t(subMatrix)[, plotRange],
        col=col[1:6], ylim=ylim, ylab=ylab, space=0, xlab=xlab, ...)

    tickPos = seq(0.5, length(plotRange)-0.5, 1)
    axis(side=1, at=tickPos, labels=FALSE, tick=TRUE, tcl=-0.3)
    for (i in 1:length(plotRange)) {
        mtext(refSeq[plotRange[i]], side=1, line=0,
            outer=FALSE, at=tickPos[i], cex=sequenceCex)
    }

    maxNumReads = ceiling(max(numReads[plotRange])/100) * 100
    scaleFac = maxPercent/maxNumReads
    lines(x=1:length(plotRange), y=numReads[plotRange]*scaleFac,
        col=col[7])
    myTicks = c(0, maxNumReads/2*scaleFac, maxNumReads*scaleFac)
    myLabels = as.character(c(0, maxNumReads/2, maxNumReads))
    axis(side=4, at=myTicks, tick=TRUE, labels=myLabels)
    mtext(side=4, line=3, "Number of reads (coverage)", las=0, col=col[7])
}


.plotVariationFrequency_file <- function(object, plotRange, ...) {

    if (file.exists(object)) {
        df = read.table(object, skip=6, header=TRUE, sep="\t", stringsAsFactors=FALSE)
        .plotVariationFrequency_data.frame(df, plotRange, ...)
    } else {
        stop(paste("File", object, "does not exist."))
    }
   
}


setMethod("plotVariationFrequency", 
    signature=signature(object="character", plotRange="numeric"),
    .plotVariationFrequency_file)
