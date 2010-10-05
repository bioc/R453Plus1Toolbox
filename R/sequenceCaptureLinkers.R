.sequenceCaptureLinkers <- function() {
  
    linkers = c(
        gSel3 = "CTCGAGAATTCTGGATCCTC",
        gSel4P = "GAGGATCCAGAATTCTCGAGTT"
    )

    linkersSet = DNAStringSet(linkers, use.names=TRUE)

    return(linkersSet)
}


.sequenceCaptureLinkers_character <- function(name) {
    
    linkers = sequenceCaptureLinkers()
    if (all(is.element(name, names(linkers)))) {
        return(linkers[is.element(names(linkers), name)])
    } else {
        stop(paste("Argument name must be on of the following:",
            paste(names(linkers), collapse=", ")))
    }
}


setMethod("sequenceCaptureLinkers", signature=signature(name="missing"),
    .sequenceCaptureLinkers)

setMethod("sequenceCaptureLinkers", signature=signature(name="character"),
    .sequenceCaptureLinkers_character)
