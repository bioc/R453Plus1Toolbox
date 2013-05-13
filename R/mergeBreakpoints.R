# private function for computing breakpoint distance
# This function assumes that both reads are from different cases.
# This can be checked for breakpoints with strand changes
# (not for breakpoints with strand chnges).
.computeBPDistance_sc <- function(chr1.1, pos1.1, chr1.2, pos1.2, 
                        chr2.1, pos2.1, chr2.2, pos2.2) {

    # inversion
    if (chr1.1 == chr1.2 & chr1.1 == chr2.1 & chr2.1 == chr2.2) {
        dist = c(abs(pos1.1 - pos2.1) + abs(pos1.2 - pos2.2),
            abs(pos1.1 - pos2.2) + abs(pos1.2 - pos2.1))
        return(min(dist))

    # not inversion
    } else if (chr1.1 == chr2.1 & chr1.2 == chr2.2) {
        return(abs(pos1.1 - pos2.1) + abs(pos1.2 - pos2.2))
    } else if (chr1.1 == chr2.2 & chr1.2 == chr2.1) {
        return(abs(pos1.1 - pos2.2) + abs(pos1.2 - pos2.1))

    # chromosomes do not fit
    } else {
        return(Inf)
    }
}


# private function for computing breakpoint distance
# This function is for breakpoints without strand changes.
.computeBPDistance_nsc <- function(strand1, chr1.5p, pos1.5p, chr1.3p, pos1.3p, 
                        strand2, chr2.5p, pos2.5p, chr2.3p, pos2.3p) {

    # data refers to the same strand
    if (strand1 == strand2) {

        if (chr1.5p == chr2.3p & chr1.3p == chr2.5p) {
            return(abs(pos1.5p - pos2.3p) + abs(pos1.3p - pos2.5p))

        } else {
            return(Inf)
        }

    # data refers to different strands
    } else {

        if (chr1.5p == chr2.5p & chr1.3p == chr2.3p) {
            return(abs(pos1.5p - pos2.5p) + abs(pos1.3p == pos2.3p))
        } else {
            return(Inf)
        }
    }
}


# private function for selecting breakpoints for merging
.markRelatedBreakpoints <- function(breakpoints, maxDist=1000) {

    #bps = seqsC1(breakpoints)
    mergeList = list()

    strandChange = sapply(commonBpsC1(breakpoints), function (x) {
        return(x[1, "strand"] != x[2, "strand"])
    })
    bpsc = names(breakpoints)[strandChange]
    bpnc = names(breakpoints)[!strandChange]


    
    # iterate through all bps with _strand changes_ and note cases
    # using this annotation:
    #
    #      Case 1        Case 2
    #
    #     3'   5'       3'   5'          
    #      |   |         |   | 
    #    - |   | +     + |   | -    
    #      | B |         | A |
    #      |   |         |   |
    #      |---|         |---| 
    #      |   |         |   | 
    #      |   |         |   |
    #    + | A | -     - | B | +
    #      |   |         |   |
    #     5'   3'       5'   3'
    #
    if (length(bpsc) > 1) {
        case = character(length(bpsc))
        names(case) = bpsc
        for (bp in bpsc) {
            ind5p =  commonBpsC1(breakpoints)[[bp]][, "localStart"] ==
                min(commonBpsC1(breakpoints)[[bp]][, "localStart"])
            if (commonBpsC1(breakpoints)[[bp]][ind5p, "strand"] == "+") {
                case[bp] = "case 1"
            } else {
                case[bp] = "case 2"
            }
        }

        # iterate breakpoints again and write merge list
        for (i in 1:(length(bpsc)-1)) {
            for (j in (i+1):length(bpsc)) {
                
                bi = bpsc[i]
                bj = bpsc[j]
                
                #different cases may be merged
                if (case[bi] != case[bj]) {
                    d = .computeBPDistance_sc(
                        chr1.1=commonBpsC1(breakpoints)[[bi]]$chr[1],
                        pos1.1=commonBpsC1(breakpoints)[[bi]]$breakpoint[1],
                        chr1.2=commonBpsC1(breakpoints)[[bi]]$chr[2],
                        pos1.2=commonBpsC1(breakpoints)[[bi]]$breakpoint[2],
                        chr2.1=commonBpsC1(breakpoints)[[bj]]$chr[1],
                        pos2.1=commonBpsC1(breakpoints)[[bj]]$breakpoint[1],
                        chr2.2=commonBpsC1(breakpoints)[[bj]]$chr[2],
                        pos2.2=commonBpsC1(breakpoints)[[bj]]$breakpoint[2])
                    
                    if (d <= maxDist) {
                        # we should merge
                        mergeBps = c(bi, bj)

                        if (any(is.element(mergeBps, unlist(mergeList)))) {
                            # we cannot merge
                            warning(paste("Want to merge ", mergeBps[1],
                                " with ", mergeBps[2], ", but one of these ",
                                "already merges with another breakpoint!",
                                sep=""))
                        } else {
                            mergeList[[length(mergeList)+1]] = mergeBps
                        }
                    }
                }
            }
        }
    }

    # Iterate through all bps with _no strand changes_ and merge bps.note
    # For same strand chimeric reads, cases 1 and 2 cannot be distinguished,
    # but when two reads are given, we can check whether they are from
    # different strands. (This is done in .computeBPDistance_nsc)
    # using this annotation:
    #
    #      Case 1        Case 2
    #
    #     3'   5'       3'   5'          
    #      |   |         |   | 
    #    + |   | -     + |   | -
    #      | B |         | A |
    #      |   |         |   |
    #      |---|         |---|
    #      |   |         |   |
    #      |   |         |   |
    #    + | A | -     + | B | -
    #      |   |         |   |
    #     5'   3'       5'   3'
    #
    if (length(bpnc) > 1) {

        # iterate breakpoints and write merge list
        for (i in 1:(length(bpnc)-1)) {
            for (j in (i+1):length(bpnc)) {

                bi = bpnc[i]
                bj = bpnc[j]
                
                ind5pi =  commonBpsC1(breakpoints)[[bi]][, "localStart"] ==
                    min(commonBpsC1(breakpoints)[[bi]][, "localStart"])
                ind5pj =  commonBpsC1(breakpoints)[[bj]][, "localStart"] ==
                    min(commonBpsC1(breakpoints)[[bj]][, "localStart"])

                d = .computeBPDistance_nsc(
                    strand1=commonBpsC1(breakpoints)[[bi]]$strand[1],
                    chr1.5p=commonBpsC1(breakpoints)[[bi]]$chr[ind5pi],
                    pos1.5p=commonBpsC1(breakpoints)[[bi]]$breakpoint[ind5pi],
                    chr1.3p=commonBpsC1(breakpoints)[[bi]]$chr[!ind5pi],
                    pos1.3p=commonBpsC1(breakpoints)[[bi]]$breakpoint[!ind5pi],
                    strand2=commonBpsC1(breakpoints)[[bj]]$strand[1],
                    chr2.5p=commonBpsC1(breakpoints)[[bj]]$chr[ind5pj],
                    pos2.5p=commonBpsC1(breakpoints)[[bj]]$breakpoint[ind5pj],
                    chr2.3p=commonBpsC1(breakpoints)[[bj]]$chr[!ind5pj],
                    pos2.3p=commonBpsC1(breakpoints)[[bj]]$breakpoint[!ind5pj])

                if (d <= maxDist) {
                    # we should merge
                    mergeBps = c(bi, bj)

                    if (any(is.element(mergeBps, unlist(mergeList)))) {
                        # we cannot merge
                        warning(paste("Want to merge ", mergeBps[1],
                            " with ", mergeBps[2], ", but one of these ",
                            "already merges with another breakpoint!",
                            sep=""))
                    } else {
                        mergeList[[length(mergeList)+1]] = mergeBps
                    }                    
                }
            }
        }
    }

    return(mergeList)
}


.mergeBreakpoints <- function(breakpoints, maxDist, mergeBPs) {

    if(length(breakpoints) == 0)
      stop("The given object contains no breakpoints")
  
    # no merging on a breakpoint object that has already been merged
    if (any(sapply(seqsC2(breakpoints), function(x) {nrow(x) != 0}))) {
        stop("The given object already contains merged breakpoints - merging aborted.")
    }

    # we need distinct breakpoint names
    if (is.null(names(breakpoints)) |
        any(is.na(names(breakpoints))) |
        any(duplicated(names(breakpoints)))) {
        stop("Breakpoints must have proper unique names.")
    }

    # compute breakpoints that should be merged
    if (missing(mergeBPs)) {
        if(missing(maxDist)) {
            maxDist = 1000
          }
        mergeBPs = .markRelatedBreakpoints(breakpoints, maxDist)
    # breakpoints to be merged are given by user
    } else {
        if (any(duplicated(unlist(mergeBPs)))) {
            stop("Breakpoints must not be merged twice!")
        }
        if (is.numeric(unlist(mergeBPs))) {
            lapply(mergeBPs, function (x) {names(breakpoints)[x]})
        }
    }

    
    # start merging
    if (length(mergeBPs) > 0) {
        for (i in 1:length(mergeBPs)) {
            k = mergeBPs[[i]][1]
            l = mergeBPs[[i]][2]

            message(paste("Merging breakpoints ", k, " and ", l, ".", sep=""))

            # seqs
            breakpoints@seqsC2[[k]] = breakpoints@seqsC1[[l]]

            # commonBps
            breakpoints@commonBpsC2[[k]] = breakpoints@commonBpsC1[[l]]

            # commonAlign (PairwiseAlignmentsSingleSubject)
            breakpoints@commonAlignC2[[k]] = breakpoints@commonAlignC1[[l]]
            
            # alignedReads (AlignedRead)
            breakpoints@alignedReadsC2[[k]] = breakpoints@alignedReadsC1[[l]]
            
            # rename  entry k
            n = names(breakpoints)
            n[which(n == k)] = paste(k, l, sep="_")
            names(breakpoints) = n

            # remove entry l
            breakpoints@seqsC1[l] = NULL
            breakpoints@seqsC2[l] = NULL
            breakpoints@commonBpsC1[l] = NULL
            breakpoints@commonBpsC2[l] = NULL
            breakpoints@commonAlignC1[l] = NULL
            breakpoints@commonAlignC2[l] = NULL
            breakpoints@alignedReadsC1[l] = NULL
            breakpoints@alignedReadsC2[l] = NULL
        }
    }

    ind = order(summary(breakpoints)$NoReadsTotal, decreasing=TRUE)
    breakpoints@seqsC1 = breakpoints@seqsC1[ind]
    breakpoints@seqsC2 =  breakpoints@seqsC2[ind]
    breakpoints@commonBpsC1 = breakpoints@commonBpsC1[ind]
    breakpoints@commonBpsC2 = breakpoints@commonBpsC2[ind]
    breakpoints@commonAlignC1 = breakpoints@commonAlignC1[ind]
    breakpoints@commonAlignC2 = breakpoints@commonAlignC2[ind]
    breakpoints@alignedReadsC1 = breakpoints@alignedReadsC1[ind]
    breakpoints@alignedReadsC2 = breakpoints@alignedReadsC2[ind]

    return(breakpoints)
}

setMethod("mergeBreakpoints",
    signature=signature(breakpoints="Breakpoints", maxDist="missing", mergeBPs="missing"),
    .mergeBreakpoints)

setMethod("mergeBreakpoints",
    signature=signature(breakpoints="Breakpoints", maxDist="numeric", mergeBPs="missing"),
    .mergeBreakpoints)

setMethod("mergeBreakpoints",
    signature=signature(breakpoints="Breakpoints", maxDist="missing", mergeBPs="list"),
    .mergeBreakpoints)
