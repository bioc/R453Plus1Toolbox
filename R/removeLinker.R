.removeLinker <- function(reads, linker,
   removeReadsWithoutLinker, minOverlap, penalty) {

    w = length(linker)
  
    if (missing(removeReadsWithoutLinker)) {
        removeReadsWithoutLinker = FALSE
    }
   
    if (missing(minOverlap)) {
        minOverlap = round(w/2)
    }
    if (missing(penalty)) {
        penalty = 2
    }
    
    readsPrefix = subseq(reads, start=1, end=pmin(width(reads), w+5))

    submat = matrix(-penalty, nrow=5, ncol=5,
        dimnames=list(c("A","C","G","T","N"), c("A","C","G","T","N")))
    diag(submat)[1:4] = 1

    pA = pairwiseAlignment(pattern=readsPrefix, subject=linker,
          type="overlap", substitutionMatrix=submat,
          gapOpening=0, gapExtension=penalty)

    ind = score(pA) >= minOverlap
    trimReads = reads
    linkerEnd = end(pA@pattern@range)+1
    linkerEnd[!ind] = 1
    trimReads = subseq(trimReads, start=linkerEnd)

    if (removeReadsWithoutLinker) {
      trimReads = trimReads[ind]
    }
    return(trimReads)
}

setMethod("removeLinker",
    signature=signature(reads="XStringSet", linker="DNAString", 
        removeReadsWithoutLinker="logical", minOverlap="numeric", penalty="numeric"),
    .removeLinker)

setMethod("removeLinker",
    signature=signature(reads="XStringSet", linker="DNAString", 
        removeReadsWithoutLinker="missing", minOverlap="missing", penalty="missing"),
    .removeLinker)
