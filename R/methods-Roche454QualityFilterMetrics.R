# Getter and setter methods
setMethod("dataTCAG", signature(object="Roche454QualityFilterMetrics"),
    function(object) {
        return(object@dataTCAG)
    }
)

setReplaceMethod("dataTCAG",
    signature=signature(object="Roche454QualityFilterMetrics",
        value="data.frame"),
    function(object, value) {
        object@dataTCAG = value
        return(object)
    }
)

setMethod("dataCATG", signature(object="Roche454QualityFilterMetrics"),
    function(object) {
        return(object@dataCATG)
    }
)

setReplaceMethod("dataCATG",
    signature=signature(object="Roche454QualityFilterMetrics",
        value="data.frame"),
    function(object, value) {
        object@dataCATG = value
        return(object)
    }
)

setMethod("dataTotal", signature(object="Roche454QualityFilterMetrics"),
    function(object) {
        return(df = cbind(object@dataTotal,
            object@dataTCAG[,c(-1,-2)] + object@dataCATG[,c(-1,-2)]))
    }
)

setReplaceMethod("dataTotal",
    signature=signature(object="Roche454QualityFilterMetrics",
        value="data.frame"),
    function(object, value) {
        object@dataTotal = value
        return(object)
    }
)


# subsequent methods are removed due to warning by R CMD check

# Subsetting
#setMethod("[", signature=signature(
#        x="Roche454QualityFilterMetrics",
#        i="ANY",
#        j="missing"),          
#    function(x, i, ..., drop) {
#        if (missing(drop)) {
#            drop = FALSE
#        }         
#        dataTotal(x) = dataTotal(x)[i, 1:3, drop=drop]
#        dataTCAG(x) = dataTCAG(x)[i, , drop=drop]
#        dataCATG(x) = dataCATG(x)[i, , drop=drop]
#        return(x)
#    }
#)


# Other convenience functions
#setMethod(show, signature(object="Roche454QualityFilterMetrics"),
#    function(object) {
#        show(dataTotal(object))
#    }
#)

#setMethod("sampleNames", signature(object="Roche454QualityFilterMetrics"),
#    function(object) {
#      return(rownames(object@dataTotal))
#    }
#)

#setMethod("combine", signature(
#        x="Roche454QualityFilterMetrics",
#        y="Roche454QualityFilterMetrics"),
#    function(x, y) {
#      if (any(is.element(sampleNames(x), sampleNames(y)))) {
#        stop(paste("Cannot combine objects: Objects contain the same regions",
#                  "from the same pico titer plate."))
#      }
#      qfm = new("Roche454QualityFilterMetrics")
#      qfm@dataTotal = rbind(x@dataTotal, y@dataTotal)
#      qfm@dataTCAG = rbind(x@dataTCAG, y@dataTCAG)
#      qfm@dataCATG = rbind(x@dataCATG, y@dataCATG)
#      return(qfm)
#    }
#)
