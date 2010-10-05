# Getter and setter methods
setMethod("dataTCAG", signature(object="Roche454BaseCallerMetrics"),
    function(object) {
        return(object@dataTCAG)
    }
)

setReplaceMethod("dataTCAG",
    signature=signature(object="Roche454BaseCallerMetrics",
        value="data.frame"),
    function(object, value) {
        object@dataTCAG = value
        return(object)
    }
)

setMethod("dataCATG", signature(object="Roche454BaseCallerMetrics"),
    function(object) {
        return(object@dataCATG)
    }
)

setReplaceMethod("dataCATG",
    signature=signature(object="Roche454BaseCallerMetrics",
        value="data.frame"),
    function(object, value) {
        object@dataCATG = value
        return(object)
    }
)

setMethod("lengthTCAG", signature(object="Roche454BaseCallerMetrics"),
    function(object) {
        return(object@lengthTCAG)
    }
)

setReplaceMethod("lengthTCAG",
    signature=signature(object="Roche454BaseCallerMetrics",
        value="matrix"),
    function(object, value) {
        object@lengthTCAG = value
        return(object)
    }
)

setMethod("scoreTCAG", signature(object="Roche454BaseCallerMetrics"),
    function(object) {
        return(object@scoreTCAG)
    }
)

setReplaceMethod("scoreTCAG",
    signature=signature(object="Roche454BaseCallerMetrics",
        value="matrix"),
    function(object, value) {
        object@scoreTCAG = value
        return(object)
    }
)

setMethod("lengthCATG", signature(object="Roche454BaseCallerMetrics"),
    function(object) {
        return(object@lengthCATG)
    }
)

setReplaceMethod("lengthCATG",
    signature=signature(object="Roche454BaseCallerMetrics",
        value="matrix"),
    function(object, value) {
        object@lengthCATG = value
        return(object)
    }
)

setMethod("scoreCATG", signature(object="Roche454BaseCallerMetrics"),
    function(object) {
        return(object@scoreCATG)
    }
)

setReplaceMethod("scoreCATG",
    signature=signature(object="Roche454BaseCallerMetrics",
        value="matrix"),
    function(object, value) {
        object@scoreTCAG = value
        return(object)
    }
)


# subsequent methods are removed due to warnings by R CMD check

# Subsetting
#setMethod("[", signature=signature(
#        x="Roche454BaseCallerMetrics",
#        i="ANY",
#        j="missing"),          
#    function(x, i, ..., drop) {
#        if (missing(drop)) {
#            drop = FALSE
#        }
#        dataTCAG(x) = dataTCAG(x)[i, , drop=drop]
#        dataCATG(x) = dataCATG(x)[i, , drop=drop]
#        lengthTCAG(x) = lengthTCAG(x)[, i, drop=drop]
#        scoreTCAG(x) = scoreTCAG(x)[, i, drop=drop]
#        lengthCATG(x) = lengthCATG(x)[, i, drop=drop]
#        scoreCATG(x) =scoreCATG(x)[, i, drop=drop]
#        return(x)
#    }
#)


# Other convenience functions
#setMethod(show, signature(object="Roche454BaseCallerMetrics"),
#    function(object) {
#        show(dataTCAG(object))
#    }
#)

#setMethod("sampleNames", signature(object="Roche454BaseCallerMetrics"),
#    function(object) {
#      return(rownames(object@dataTCAG))
#    }
#)

#setMethod("combine", signature(
#        x="Roche454BaseCallerMetrics",
#        y="Roche454BaseCallerMetrics"),
#    function(x, y) {
#      if (any(is.element(sampleNames(x), sampleNames(y)))) {
#        stop(paste("Cannot combine objects: Objects contain the same regions",
#                   "from the same pico titer plate."))
#      }
#      bcm = new("Roche454BaseCallerMetrics")
#      bcm@dataTCAG = rbind(x@dataTCAG, y@dataTCAG)
#      bcm@dataCATG = rbind(x@dataCATG, y@dataCATG)
#
#      ncolX = ncol(x@lengthTCAG)
#      ncolY = ncol(y@lengthTCAG)
#      z = matrix(0, nrow=max(nrow(x@lengthTCAG), nrow(y@lengthTCAG)),
#          ncol=ncolX+ncolY)
#      z[1:nrow(x@lengthTCAG), 1:ncolX] = x@lengthTCAG
#      z[1:nrow(y@lengthTCAG), (ncolX+1):(ncolX+ncolY)] = y@lengthTCAG
#      rownames(z) = c(rownames(x@lengthTCAG), rownames(y@lengthTCAG))
#      bcm@lengthTCAG = z
#
#      ncolX = ncol(x@scoreTCAG)
#      ncolY = ncol(y@scoreTCAG)
#      z = matrix(0, nrow=max(nrow(x@scoreTCAG), nrow(y@scoreTCAG)),
#          ncol=ncolX+ncolY)
#      z[1:nrow(x@scoreTCAG), 1:ncolX] = x@scoreTCAG
#      z[1:nrow(y@scoreTCAG), (ncolX+1):(ncolX+ncolY)] = y@scoreTCAG
#      rownames(z) = c(rownames(x@scoreTCAG), rownames(y@scoreTCAG))
#      bcm@scoreTCAG = z
#
#      ncolX = ncol(x@lengthCATG)
#      ncolY = ncol(y@lengthCATG)
#      z = matrix(0, nrow=max(nrow(x@lengthCATG), nrow(y@lengthCATG)),
#          ncol=ncolX+ncolY)
#      z[1:nrow(x@lengthCATG), 1:ncolX] = x@lengthCATG
#      z[1:nrow(y@lengthCATG), (ncolX+1):(ncolX+ncolY)] = y@lengthCATG
#      rownames(z) = c(rownames(x@lengthCATG), rownames(y@lengthCATG))
#      bcm@lengthCATG = z
#
#      ncolX = ncol(x@scoreCATG)
#      ncolY = ncol(y@scoreCATG)
#      z = matrix(0, nrow=max(nrow(x@scoreCATG), nrow(y@scoreCATG)),
#          ncol=ncolX+ncolY)
#      z[1:nrow(x@scoreCATG), 1:ncolX] = x@scoreCATG
#      z[1:nrow(y@scoreCATG), (ncolX+1):(ncolX+ncolY)] = y@scoreCATG
#      rownames(z) = c(rownames(x@scoreCATG), rownames(y@scoreCATG))
#      bcm@scoreCATG = z
#
#      return(bcm)
#    }
#)
