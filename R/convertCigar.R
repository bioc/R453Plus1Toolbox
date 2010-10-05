# These are temporary methods, that are likely to be replaced by
# methods from the Rsamtools package in near future.


# takes a character vector with CIGAR strings, returns a list
extendedCIGARToList <- function(cigars) {
    indList = gregexpr("[A-Z]", cigars)
    cigarList = list()

    for (l in 1:length(indList)) {
        ind = c(0, indList[[l]])
        s = cigars[[l]]
        m = integer()
        for (i in 2:length(ind)) {
            m[i-1] = as.integer(substr(s, ind[i-1]+1, ind[i]-1))
            names(m)[i-1] = substr(s, ind[i], ind[i])
        }
        cigarList[[l]] = m
    }
    names(cigarList) = names(cigars)
    return(cigarList)
}

# takes a list with converted CIGAR string (see above) and returns a character vector
listToExtendedCIGAR <- function(cigarList){

    cigars = vector(mode="character")

    for (l in 1:length(cigarList))
	 cigars[l] = paste(paste(cigarList[[l]], names(cigarList[[l]]), sep=""), collapse="")

    names(cigars) = names(cigarList)
    return(cigars)
}
