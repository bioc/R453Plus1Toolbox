.plotAmpliconCoverage <- function(avaSet, type="amplicon", bothDirections=TRUE,
    cex.names, cex.axis, las, col, ...) {
  
    args = list(...)
   
##########################################
# Plot number of reads for each amplicon #
##########################################
    if (type == "amplicon") {

        if (ncol(avaSet) == 1) {

            # default values for main and ylab
            main = sampleNames(avaSet)
            ylab = "Number of reads (coverage)"

                 
            data = cbind(assayDataAmp(avaSet)$forwCount,
                assayDataAmp(avaSet)$revCount)
            data = data[sort(rownames(data)), ]
            data = t(data)

            maxCov = max(data[, 1] + data[, 2])
            maxCov = maxCov + (0.2 * maxCov)    ## add some extra space for the legend

            par(oma=c(4,1,0,0))
            if (is.element("main", names(args)) & is.element("ylab", names(args))) {
                barplot(data, ylim=c(0, maxCov),
                    cex.names=cex.names, names.arg=colnames(data), las=las, col=col, ...)
            } else if (is.element("main", names(args)) & !is.element("ylab", names(args))) {
                barplot(data, ylab=ylab, ylim=c(0, maxCov),
                    cex.names=cex.names, names.arg=colnames(data), las=las, col=col, ...)
            } else if (!is.element("main", names(args)) & is.element("ylab", names(args))) {
                barplot(data, main=main, ylim=c(0, maxCov),
                    cex.names=cex.names, names.arg=colnames(data), las=las, col=col, ...)
            } else {
                barplot(data, ylab=ylab, main=main, ylim=c(0, maxCov),
                    cex.names=cex.names, names.arg=colnames(data), las=las, col=col, ...)
            }
            legend(x="topright", legend=c("forward reads", "reverse reads"), fill=col)       

            
        } else {

            # default values for main and ylab
            main = "Amplicon coverage"
            ylab = "Number of reads (coverage)"

            numAmps = nrow(assayDataAmp(avaSet)$forwCount)
            ampNames = sort(rownames(assayDataAmp(avaSet)$forwCount))

            if (bothDirections) {
                              
                data = matrix(NA, nrow=ncol(avaSet), ncol=2*numAmps)
                colnames(data) = rep("", numAmps*2)
                data[, seq(1, 2*numAmps, 2)] =
                    t(assayDataAmp(avaSet)$forwCount[ampNames, ])
                data[, seq(2, 2*numAmps, 2)] =
                    t(assayDataAmp(avaSet)$revCount[ampNames, ])
                colnames(data)[seq(1, 2*numAmps, 2)] = ampNames
                colnames(data)[seq(2, 2*numAmps, 2)] = ampNames

                maxCov = max(data)
                maxCov = maxCov + (0.2 * maxCov)    ## add some extra space for the legend

                par(oma=c(4,1,0,0))
                if (is.element("main", names(args)) & is.element("ylab", names(args))) {
                    boxplot(data, las=las, ylim=c(0, maxCov),
                        cex.axis=cex.axis, col=col, ...)
                } else if (is.element("main", names(args)) & !is.element("ylab", names(args))) {
                    boxplot(data, las=las, ylim=c(0, maxCov), ylab=ylab,
                        cex.axis=cex.axis, col=col, ...)
                } else if (!is.element("main", names(args)) & is.element("ylab", names(args))) {
                    boxplot(data, las=las, main=main, ylim=c(0, maxCov),
                        cex.axis=cex.axis, col=col, ...)
                } else {
                    boxplot(data, las=las, main=main, ylim=c(0, maxCov), ylab=ylab,
                        cex.axis=cex.axis, col=col, ...)
                }
                legend("topright", legend=c("forward reads", "reverse reads"), fill=col)
                
            } else {
                
                data = t(assayDataAmp(avaSet)$forwCount[ampNames, ]) +
                    t(assayDataAmp(avaSet)$revCount[ampNames, ])
                colnames(data) = ampNames

                par(oma=c(4,1,0,0))
                if (is.element("main", names(args)) & is.element("ylab", names(args))) {
                    boxplot(data, las=las,
                    cex.axis=cex.axis, col=col[1], ...)
                } else if (is.element("main", names(args)) & !is.element("ylab", names(args))) {
                    boxplot(data, las=las, ylab=ylab,
                    cex.axis=cex.axis, col=col[1], ...)
                } else if (!is.element("main", names(args)) & is.element("ylab", names(args))) {
                    boxplot(data, las=las, main=main,
                    cex.axis=cex.axis, col=col[1], ...)
                } else {
                    boxplot(data, las=las, main=main, ylab=ylab,
                    cex.axis=cex.axis, col=col[1], ...)
                }
            } 
        }


#####################################
# Plot number of reads for each MID #
#####################################
    } else if (type == "mid") {

        #default values for main and ylab
        main = "Number of reads by MID"
        ylab = "Number of reads (per sample)"
      
        if (any(!is.na(pData(avaSet)$MID2)) &
            any(pData(avaSet)$MID1 != pData(avaSet)$MID2)) {
            stop("MID1 and MID2 must be equal.")
        }
        ind = grep(",", pData(avaSet)$MID1, fixed=TRUE)
        if (length(ind) != 0) {
            avaSet = avaSet[, -ind]
            warning(paste(length(ind), "samples were excluded due to mutiple MIDs."))
        }

        if (bothDirections) {
          
            midFac = factor(c(paste(as.character(pData(avaSet)$MID1), "_f", sep=""),
                paste(as.character(pData(avaSet)$MID1), "_r", sep="")))
            data = c(apply(assayDataAmp(avaSet)$forwCount, 2, sum),
                apply(assayDataAmp(avaSet)$revCount, 2, sum))

            maxCov = max(data)
            maxCov = maxCov + (0.2 * maxCov)    ## add some extra space for the legend

            par(oma=c(2,1,0,0))
            if (is.element("main", names(args)) & is.element("ylab", names(args))) {
                boxplot(data ~ midFac, las=las, col=col, ylim=c(0, maxCov),
                    names=substring(levels(midFac), first=1, last=nchar(levels(midFac))-2),
                    cex.axis=cex.axis, ...)
            } else if (is.element("main", names(args)) & !is.element("ylab", names(args))) {
                boxplot(data ~ midFac, las=las, col=col, ylim=c(0, maxCov), ylab=ylab,
                    names=substring(levels(midFac), first=1, last=nchar(levels(midFac))-2),
                    cex.axis=cex.axis, ...)
            } else if (!is.element("main", names(args)) & is.element("ylab", names(args))) {
                boxplot(data ~ midFac, las=las, col=col, main=main, ylim=c(0, maxCov),
                    names=substring(levels(midFac), first=1, last=nchar(levels(midFac))-2),
                    cex.axis=cex.axis, ...)
            } else {
                boxplot(data ~ midFac, las=las, col=col, main=main, ylim=c(0, maxCov), ylab=ylab,
                    names=substring(levels(midFac), first=1, last=nchar(levels(midFac))-2),
                    cex.axis=cex.axis, ...)
            }
            legend("topright", legend=c("forward reads", "reverse reads"), fill=col)

        } else {
          
            midFac = factor(as.character(pData(avaSet)$MID1))
            data = apply(assayDataAmp(avaSet)$forwCount, 2, sum) +
                apply(assayDataAmp(avaSet)$revCount, 2, sum)

            par(oma=c(2,1,0,0))
            if (is.element("main", names(args)) & is.element("ylab", names(args))) {
                boxplot(data ~ midFac, las=las, col=col[1], 
                    cex.axis=cex.axis, ...)
            } else if (is.element("main", names(args)) & !is.element("ylab", names(args))) {
                boxplot(data ~ midFac, las=las, col=col[1], ylab=ylab,
                    cex.axis=cex.axis, ...)
            } else if (!is.element("main", names(args)) & is.element("ylab", names(args))) {
                boxplot(data ~ midFac, las=las, col=col[1], main=main, 
                    cex.axis=cex.axis, ...)
            } else {
                boxplot(data ~ midFac, las=las, col=col[1], main=main, ylab=ylab,
                    cex.axis=cex.axis, ...)
            }

        }


##############################################
# Plot Number of reads for each PTP / region #
##############################################
    } else if (type == "ptp") {

        # default values for main and ylab
        main = "Number of reads by pico titer plate / region"
        ylab = "Number of reads (per sample)"
      
        ind = union(grep(",", pData(avaSet)$PTP_AccNum, fixed=TRUE),
            grep(",", pData(avaSet)$Lane, fixed=TRUE))
        if (length(ind) != 0) {
            avaSet = avaSet[, -ind]
            warning(paste(length(ind), "samples were excluded due to mutiple PTPs/Lanes."))
        }

        if (bothDirections) {
            ptpFac = paste(as.character(pData(avaSet)$PTP_AccNum),
                as.character(pData(avaSet)$Lane), sep="")
            ptpFac = factor(c(paste(ptpFac, "_f", sep=""), paste(ptpFac, "_r", sep="")))
            data = c(apply(assayDataAmp(avaSet)$forwCount, 2, sum),
                apply(assayDataAmp(avaSet)$revCount, 2, sum))

            maxCov = max(data)
            maxCov = maxCov + (0.2 * maxCov)    ## add some extra space for the legend

            par(oma=c(2,1,0,0))

            if (is.element("main", names(args)) & is.element("ylab", names(args))) {
                boxplot(data ~ ptpFac, las=las, col=col, ylim=c(0, maxCov),
                    names=substring(levels(ptpFac), first=1, last=nchar(levels(ptpFac))-2),
                    cex.axis=cex.axis, ...)
            } else if (is.element("main", names(args)) & !is.element("ylab", names(args))) {
                boxplot(data ~ ptpFac, las=las, col=col, ylim=c(0, maxCov), ylab=ylab,
                    names=substring(levels(ptpFac), first=1, last=nchar(levels(ptpFac))-2),
                    cex.axis=cex.axis, ...)
            } else if (!is.element("main", names(args)) & is.element("ylab", names(args))) {
                boxplot(data ~ ptpFac, las=las, col=col, main=main, ylim=c(0, maxCov),
                    names=substring(levels(ptpFac), first=1, last=nchar(levels(ptpFac))-2),
                    cex.axis=cex.axis, ...)
            } else {
                boxplot(data ~ ptpFac, las=las, col=col, main=main, ylim=c(0, maxCov), ylab=ylab,
                    names=substring(levels(ptpFac), first=1, last=nchar(levels(ptpFac))-2),
                    cex.axis=cex.axis, ...)
            }
            legend("topright", legend=c("forward reads", "reverse reads"), fill=col)

        } else {
          
            ptpFac = factor(paste(as.character(pData(avaSet)$PTP_AccNum),
                as.character(pData(avaSet)$Lane), sep=""))
            data = apply(assayDataAmp(avaSet)$forwCount, 2, sum) +
                apply(assayDataAmp(avaSet)$revCount, 2, sum)

            par(oma=c(2,1,0,0))
            if (is.element("main", names(args)) & is.element("ylab", names(args))) {
                boxplot(data ~ ptpFac, las=las, col=col[1], 
                    cex.axis=cex.axis, ...)
            } else if (is.element("main", names(args)) & !is.element("ylab", names(args))) {
                boxplot(data ~ ptpFac, las=las, col=col[1], ylab=ylab,
                    cex.axis=cex.axis, ...)
            } else if (!is.element("main", names(args)) & is.element("ylab", names(args))) {
                boxplot(data ~ ptpFac, las=las, col=col[1], main=main,
                    cex.axis=cex.axis, ...)
            } else {
                boxplot(data ~ ptpFac, las=las, col=col[1], main=main, ylab=ylab,
                    cex.axis=cex.axis, ...)
            }
        }


    } else {
        print("Argument type must be amplicon, mid or ptp.")
    }

}


setMethod("plotAmpliconCoverage", 
    signature=signature(avaSet="AVASet", type="character", bothDirections="logical"),
    .plotAmpliconCoverage)


setMethod("plotAmpliconCoverage",
    signature=signature(avaSet="AVASet", type="character", bothDirections="missing"),
    .plotAmpliconCoverage)


setMethod("plotAmpliconCoverage",
    signature=signature(avaSet="AVASet", type="missing", bothDirections="logical"),
    .plotAmpliconCoverage)
    


setMethod("plotAmpliconCoverage",
    signature=signature(avaSet="AVASet", type="missing", bothDirections="missing"),
    .plotAmpliconCoverage)
