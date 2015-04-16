.demultiplexReads <- function(reads, mids, numMismatches, trim) {

    if (missing(numMismatches)) {
        numMismatches = 2
    }
    if (missing(trim)) {
        trim = TRUE
    }
    
    w = width(mids)[1]
    if (!all(width(mids) == w)) {
        stop("All MIDs must have the same length.")
    }
    
    readsPrefix = subseq(reads, start=1, end=pmin(width(reads), w+5))

    submat = matrix(-1, nrow=5, ncol=5,
        dimnames=list(c("A","C","G","T","N"), c("A","C","G","T","N")))
    diag(submat)[1:4] = 0

    scores = matrix(NA, nrow=length(readsPrefix), ncol=length(mids))
    patternRanges = list()
    for (i in 1:length(mids)) {
        pA = pairwiseAlignment(pattern=readsPrefix, subject=mids[[i]],
            type="local-global", substitutionMatrix=submat,
            gapOpening=0, gapExtension=1)
        scores[,i] = score(pA)
        patternRanges[[i]] = pA@pattern@range
    }

    # find unique maxima >= numMismatches
    sample = apply(scores, 1,
        function(x) {
            m = max(x)
            if (sum(x == m) == 1 & m >= -numMismatches) {
                return(which(x == m))
            } else {
                return(NA)
            }
        }
    )

    sampleSeqs = list()
    for (i in 1:length(mids)) {
        ind = sample == i
        ind[is.na(ind)] = FALSE
        if (trim & any(ind)) {
           sampleSeqs[[i]] = subseq(reads[ind],
               start=end(patternRanges[[i]][ind])+1)
        } else {
           sampleSeqs[[i]] = reads[ind]
        }    
    }
    names(sampleSeqs) = names(mids)

    return(sampleSeqs)
}

setMethod("demultiplexReads",
    signature=signature(reads="XStringSet", mids="XStringSet", numMismatches="numeric", trim="logical"),
    .demultiplexReads)

setMethod("demultiplexReads",
    signature=signature(reads="XStringSet", mids="XStringSet", numMismatches="numeric", trim="missing"),
    .demultiplexReads)

setMethod("demultiplexReads",
    signature=signature(reads="XStringSet", mids="XStringSet", numMismatches="missing", trim="logical"),
    .demultiplexReads)

setMethod("demultiplexReads",
    signature=signature(reads="XStringSet", mids="XStringSet", numMismatches="missing", trim="missing"),
    .demultiplexReads)

