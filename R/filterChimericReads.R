.filterChimericReads <- function(alnReads, targetRegion, linkerSeq,
                                 minDist, dupReadDist) {

    if (missing(minDist)) {
        minDist = 1000
    }
    if (missing(dupReadDist)) {
        dupReadDist = 1
    }

    # check sam/bam input
    reads = lapply(names(alnReads[[1]]), function(elt) {
        if (is.factor(alnReads[[1]][[elt]])) {
            factor(do.call(c, unname(lapply(alnReads, function (x) {
                as.character(x[[elt]])
            }))))
                                            
        } else
            do.call(c, unname(lapply(alnReads, "[[", elt)))
        }
    )
    names(reads) = names(alnReads[[1]])
    reqFields = c("qname", "rname", "pos", "cigar")
    if (!all(is.element(reqFields, names(reads)))) {
        stop(paste("The following SAM fields are necessary:", reqFields))
    }
    reads = do.call("DataFrame", reads)


    # initiate data frame for logging info
    log = data.frame(AlignedReads=length(unique(reads$qname)),
        stringsAsFactors=FALSE)

    
    # 1: extract reads with exactly two local alignments
    message("Extrating chimeric reads with exactly two local alignments.")
    rTable = table(reads$qname)
    crNames = names(rTable[rTable == 2])
    reads = reads[is.element(reads$qname, crNames),]
    log$ChimericReads = sum(rTable > 1)
    log$TwoLocalAlignments = sum(rTable == 2)

    
    # intermediate step : add these four columns to the data frame
    # start      - start position of the alignment in reference coordinates
    #              (equal to the SAM field pos)
    # end        - end position of the alignment in reference coordinates
    #              (corrected for insertions in the query read)
    # localStart - start position of the alignment in query seq. coordinates
    # localEnd   - end position of the alignment in query seq. coordinates
    #              (corrected for deletions in the query read)
    # All fields refer to the positive strand as usual in SAM.
    cigars = extendedCIGARToList(reads$cigar)
    dels = sapply(cigars, function(x) sum(x[names(x)=="D"]))
    ins = sapply(cigars, function(x) sum(x[names(x)=="I"]))
    matches = sapply(cigars, function(x) sum(x[names(x)=="M"]))
    reads$start = reads$pos
    reads$end = as.integer(reads$start + matches + dels - 1)
    reads$localStart = as.integer(sapply(cigars, FUN=
        function(x) {
            ind = min(which(names(x) == "M"))
            x = c(0, x)
            return(sum(x[1:ind]) + 1)
        }
    ))
    reads$localEnd = as.integer(sapply(cigars, FUN=
        function(x) {
            ind = max(which(names(x) == "M"))
            return(sum(x[1:ind]))
        }
    ) - dels)
    

    # 2: remove all reads that do not align to the chip's target region
    if (!missing(targetRegion)) {
        message("Removing reads that do not align to the target region.")
        readRanges = GRanges(
            IRanges(start=reads$start, end=reads$end),
            seqnames=reads$rname, name=reads$qname)
        ind = suppressWarnings(overlapsAny(readRanges, GRanges(targetRegion)))
        targetReads <- readRanges$name[ind]
        reads = reads[is.element(reads$qname, targetReads),]
        log$TargetRegion = length(unique(reads$qname))
    } else {
        targetRegion = "missing" # only for parameter slot in the returned object
    }
    reads$seq = compact(reads$seq)


    # 3: remove all reads that have a linker sequence in the middle
    if (!missing(linkerSeq)) {
        message("Removing reads with a linker sequence in the middle.")
        dupInd = duplicated(reads$qname)
        seqs = reads$seq[!dupInd]
        names(seqs) = reads$qname[!dupInd]
        readsLinker = character()

        f = vmatchPattern(linkerSeq, seqs, max.mismatch=4)
        indF = elementNROWS(f) > 0
        seqLength = width(seqs)[indF]
        for (i in setdiff(0:sum(indF), 0)) {
            if (! (any(startIndex(f)[indF][[i]] < 10) |
                any(seqLength[i] - endIndex(f)[indF][[i]] < 10))) {
                readsLinker = c(readsLinker, names(seqs)[indF][i])
            }
        }
      
        r = vmatchPattern(linkerSeq, reverseComplement(seqs), max.mismatch=4)
        indR = elementNROWS(r) > 0
        seqLength = width(seqs)[indR]
        for (i in setdiff(0:sum(indR), 0)) {
            if (! (any(startIndex(r)[indR][[i]] < 10) |
                any(seqLength[i] - endIndex(r)[indR][[i]] < 10))) {
                readsLinker = c(readsLinker, names(seqs)[indR][i])
            }
        }
    
        ind = !is.element(reads$qname, readsLinker)
        reads = reads[ind, ]
        log$NoLinker = length(unique(reads$qname))
    } else {
        linkerSeq = "missing" # only for parameter slot in the returned object
    }


    # 4: Remove local alignments that are too close together (minDist)
    message("Removing reads with close local alignments (perhaps small indels).")
    cList = list()
    readNames = unique(reads$qname)
    for (n in readNames) {
        ind = reads$qname == n
        if (sum(ind) == 2) {
            if (reads[ind, "rname"][1] == reads[ind, "rname"][2]) {
                # Reads from same chr. must not overlap in reference seq.
                # and must insure a minimal distance minDist
                ir = IRanges(start=reads[ind, "start"],
                    end=reads[ind, "end"])
                if (countOverlaps(ir[1], ir[2]) == 0  &
                    min(abs(start(ir[1]) - end(ir[2])),
                        abs(end(ir[1]) - start(ir[2]))) >= minDist) {
                    cList[[n]] = reads[ind, ]
                    cList[[n]]$seq = compact(cList[[n]]$seq)
                }
            } else {
                cList[[n]] = reads[ind, ]
                cList[[n]]$seq = compact(cList[[n]]$seq)
            }
        }
    }
    log$MinimumDistance = length(cList)


    # intermediate step - extract Breakpoints:
    # Each chimeric read with two local alignments has 4 points where an
    # alignment starts/ends. The two points in the middle are the breakpoint.
    # Before we can compute the breakpoint, we need the alignment positions
    # with respect to the reads:
    # Positions are given for the sequences stored in the sam file. Seqs that
    # were aligned to the minus strand of the reference are stored as reverse
    # complement. Here, we compute start (end) and localStart (localEnd)
    # coordinates by defining the 5 prime end as start point. So, two local
    # alignments of the same read are defined in the same coordinates even if
    # one part maps to the + and the other to the - strand.
    cList = lapply(cList,
        function (df) {
            ind = df$strand == "-"
            df$localStart5 = df$localStart
            df$localEnd5 = df$localEnd
            df$start5 = df$start
            df$end5 = df$end
            if (sum(ind) > 0) {
                df$localEnd5[ind] = (df$qwidth[ind] + 1) - df$localStart[ind]
                df$localStart5[ind] = (df$qwidth[ind] + 1) - df$localEnd[ind]
                df$start5[ind] = df$end[ind]
                df$end5[ind] = df$start[ind]
            }

            bp1 = which(df$localStart5 == max(df$localStart5))[1]
            bp2 = which(df$localEnd5 == min(df$localEnd5))[1]
            
            if (bp1 == 1 & bp2 == 2) {
                df$breakpoint = c(df$start5[bp1], df$end5[bp2])
            } else if (bp1 == 2 & bp2 == 1) {
                df$breakpoint = c(df$end5[bp2], df$start5[bp1])
            } else {
                df$breakpoint = NA
            }
            return(df)
        }
    )


    # 5: Remove duplicated reads,
    # (dup. reads := reads with same chr., same strand and same 5' start pos.
    message("Removing duplicated reads (perhaps duplication due to PCR).")
    checkDF = data.frame(
        chr = sapply(cList, function(x) {
            ind = which(x$localStart5 == min(x$localStart5))[1]
            return(as.character(x$rname[ind]))
        }),
        strand = sapply(cList, function(x) {
            ind = which(x$localStart5 == min(x$localStart5))[1]
            return(as.character(x$strand[ind]))
        }),
        pos = sapply(cList, function(x) {
            ind = which(x$localStart5 == min(x$localStart5))[1]
            return(x$start5[ind])
        }),
        qwidth = sapply(cList, function(x) {x$qwidth[1]})
    )
    row.names(checkDF) = names(cList)
    dupNames = character()
    sCheckDF = split(checkDF, checkDF$chr)
    for (i in 1:length(sCheckDF)) {
        df = split(sCheckDF[[i]], sCheckDF[[i]]$strand)
        for (j in 1:length(df)) {
            if (nrow(df[[j]]) > 0) {
                for (k in 1:nrow(df[[j]])) {
                    ind = abs(df[[j]]$pos[k] - df[[j]]$pos[-k]) <= dupReadDist
                    allNames = sort(
                        c(row.names(df[[j]][k,]), row.names(df[[j]][-k,])[ind]))                    
                    longestRead = which(df[[j]][allNames, "qwidth"] ==
                        max(df[[j]][allNames, "qwidth"]))
                    dupNames = union(dupNames, allNames[-longestRead])
                }
            }
        }
    }
    cList = cList[setdiff(names(cList), dupNames)]
    log$Unique5PrimeStart = length(cList)


    # return list with filtered chimeric reads
    filteredChimReads = names(cList)
    for (i in 1:length(alnReads)) {
        ind = is.element(alnReads[[i]]$qname, filteredChimReads)
        for (j in 1:length(alnReads[[i]])) {
            alnReads[[i]][[j]] = alnReads[[i]][[j]][ind]
        }
    }
    alnReads[["log"]] = log
    return(alnReads)
}


setMethod("filterChimericReads",
    signature=signature(alnReads="list", targetRegion="IntegerRangesList", linkerSeq="DNAString",
        minDist="numeric", dupReadDist="numeric"),
    .filterChimericReads)

setMethod("filterChimericReads",
    signature=signature(alnReads="list", targetRegion="IntegerRangesList", linkerSeq="DNAString",
        minDist="missing", dupReadDist="missing"),
    .filterChimericReads)

setMethod("filterChimericReads",
    signature=signature(alnReads="list", targetRegion="IntegerRangesList", linkerSeq="missing",
        minDist="numeric", dupReadDist="numeric"),
    .filterChimericReads)

setMethod("filterChimericReads",
    signature=signature(alnReads="list", targetRegion="IntegerRangesList", linkerSeq="missing",
        minDist="missing", dupReadDist="missing"),
    .filterChimericReads)

setMethod("filterChimericReads",
    signature=signature(alnReads="list", targetRegion="missing", linkerSeq="DNAString",
        minDist="numeric", dupReadDist="numeric"),
    .filterChimericReads)

setMethod("filterChimericReads",
    signature=signature(alnReads="list", targetRegion="missing", linkerSeq="DNAString",
        minDist="missing", dupReadDist="missing"),
    .filterChimericReads)

setMethod("filterChimericReads",
    signature=signature(alnReads="list", targetRegion="missing", linkerSeq="missing",
        minDist="numeric", dupReadDist="numeric"),
    .filterChimericReads)

setMethod("filterChimericReads",
    signature=signature(alnReads="list", targetRegion="missing", linkerSeq="missing",
        minDist="missing", dupReadDist="missing"),
    .filterChimericReads)

