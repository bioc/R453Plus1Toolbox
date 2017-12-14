.coverageOnTarget <- function(alnReads, targetRegion) {

    result = list()
    for (i in 1:length(alnReads)) {
        rd = GRanges(IRanges(start=alnReads[[i]]$pos,
                             end=alnReads[[i]]$pos + alnReads[[i]]$qwidth - 1),
                             seqnames=alnReads[[i]]$rname)
        cov = coverage(rd)

        chrs = as.character(intersect(names(cov), seqnames(targetRegion)))
        result[[i]] = numeric()
        for (chr in chrs) {
            ind = end(targetRegion[seqnames(targetRegion) == chr]) > sum(width(cov[[chr]]))
            end(targetRegion[seqnames(targetRegion) == chr])[ind] = sum(width(cov[[chr]]))
            tmpTC = unlist(aggregate(cov[[chr]], targetRegion[seqnames(targetRegion) == chr], FUN=as.numeric))
            names(tmpTC) = rep(chr, length(tmpTC))
            result[[i]] = as.integer(c(result[[i]], tmpTC))
        }
    }

    return(result)
}


setMethod("coverageOnTarget", 
    signature=signature(alnReads="list", targetRegion="GRanges"),
    .coverageOnTarget)



# targetRegion muss Klasse "RangedData" sein
# readStatus ab<C3><A4>ndern in GSMSet
#
#targetCoverage <- function(gsmSet, targetRegion) {
#
#    first = TRUE
#    readStatus = getReadStatus(gsmSet)
#    
#    for (i in 1:length(readStatus)) {
#        rStat = readStatus[[i]]
#        rStat = rStat[rStat$Status == "Full", ]
#
#        arRanges = RangedData(
#            ranges=IRanges(start=rStat$Start, end=rStat$Stop),
#            space=rStat$Accno.1)
#
#        chrs = intersect(names(arRanges), names(targetRegion))
#        tRanges = targetRegion[chrs]
#        rRanges = arRanges[chrs]
#
#        cov = countOverlaps(query=rRanges, subject=tRanges)
#        sov = subsetByOverlaps(query=rRanges, subject=tRanges)
#
#        targetCov = list()
#        for (i in 1:length(chrs)) {
#            targetCov[[i]] = coverage(rRanges[chrs[i]])[[chrs[i]]]
#            targetCov[[i]] = targetCov[[i]][ranges(tRanges[chrs[i]])[[chrs[i]]]]
#        }
#        names(targetCov) = chrs
#
#        df = data.frame(
#            NoMappedReads=nrow(rStat),
#            NoTargetReads=sum(sapply(cov, sum)),
#            NoMappedBases=sum(width(arRanges)),
#            NoTargetBases=sum(sapply(targetCov, sum)),
#            MeanTargetCoverage=mean(unlist(lapply(targetCov, as.integer))),
#            MedianTargetCoverage=median(unlist(lapply(targetCov, as.integer))))
#
#        if (first) {
#            summaryStats = df
#            first = FALSE
#        } else {
#            summaryStats = rbind(summaryStats, df)
#        }
#    }
#
#    rownames(summaryStats) = names(readStatus)
#    return(summaryStats)
#}
