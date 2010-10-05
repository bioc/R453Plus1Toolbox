.genomeSequencerMIDs <- function() {
  
    mids = c(
        MID1 = "ACGAGTGCGT",
        MID2 = "ACGCTCGACA",
        MID3 = "AGACGCACTC",
        MID4 = "AGCACTGTAG",
        MID5 = "ATCAGACACG",
        MID6 = "ATATCGCGAG",
        MID7 = "CGTGTCTCTA",
        MID8 = "CTCGCGTGTC",
        MID9 = "TAGTATCAGC",
        MID10 = "TCTCTATGCG",
        MID11 = "TGATACGTCT",
        MID12 = "TACTGAGCTA",
        MID13 = "CATAGTAGTG",
        MID14 = "CGAGAGATAC"
    )

    midsSet = DNAStringSet(mids, use.names=TRUE)

    return(midsSet)
}


.genomeSequencerMIDs_character <- function(mid) {
    
    mids = genomeSequencerMIDs()
    if (all(is.element(mid, names(mids)))) {
        return(mids[is.element(names(mids), mid)])
    } else {
        stop(paste("Argument mid must be on of the following:",
            paste(names(mids), collapse=", ")))
    }
}


setMethod("genomeSequencerMIDs", signature=signature(mid="missing"),
    .genomeSequencerMIDs)

setMethod("genomeSequencerMIDs", signature=signature(mid="character"),
    .genomeSequencerMIDs_character)
