### Getter & Setter ###
setMethod("readname", signature(object="SFFRead"),
  function(object) {
    return(object@readname)
  }
)
setReplaceMethod("readname", signature(object="SFFRead", value="character"),
  function(object, value) {
    object@readname = value
    return(object)
  }
)

setMethod("bases", signature(object="SFFRead"),
  function(object) {
    return(object@bases)
  }
)
setReplaceMethod("bases", signature(object="SFFRead", value="DNAString"),
  function(object, value) {
    object@bases = value
    return(object)
  }
)

setMethod("clipQualityLeft", signature(object="SFFRead"),
  function(object) {
    return(object@clipQualityLeft)
  }
)
setReplaceMethod("clipQualityLeft", signature(object="SFFRead", value="numeric"),
  function(object, value) {
    object@clipQualityLeft = value
    return(object)
  }
)

setMethod("clipQualityRight", signature(object="SFFRead"),
  function(object) {
    return(object@clipQualityRight)
  }
)
setReplaceMethod("clipQualityRight", signature(object="SFFRead", value="numeric"),
  function(object, value) {
    object@clipQualityRight = value
    return(object)
  }
)

setMethod("clipAdapterLeft", signature(object="SFFRead"),
  function(object) {
    return(object@clipAdapterLeft)
  }
)
setReplaceMethod("clipAdapterLeft", signature(object="SFFRead", value="numeric"),
  function(object, value) {
    object@clipAdapterLeft = value
    return(object)
  }
)

setMethod("clipAdapterRight", signature(object="SFFRead"),
  function(object) {
    return(object@clipAdapterRight)
  }
)
setReplaceMethod("clipAdapterRight", signature(object="SFFRead", value="numeric"),
  function(object, value) {
    object@clipAdapterRight = value
    return(object)
  }
)

setMethod("flowgramValues", signature(object="SFFRead"),
  function(object) {
    return(object@flowgramValues)
  }
)
setReplaceMethod("flowgramValues", signature(object="SFFRead", value="numeric"),
  function(object, value) {
    object@flowgramValues = value
    return(object)
  }
)

setMethod("flowIndexes", signature(object="SFFRead"),
  function(object) {
    return(object@flowIndexes)
  }
)
setReplaceMethod("flowIndexes", signature(object="SFFRead", value="numeric"),
  function(object, value) {
    object@flowIndexes = value
    return(object)
  }
)

#setMethod("qualityScores", signature(object="SFFRead"),
#  function(object) {
#    return(object@qualityScores)
#  }
#)
#setReplaceMethod("qualityScores", signature(object="SFFRead", value="numeric"),
#  function(object, value) {
#    object@qualityScores = value
#    return(object)
#  }
#)

setMethod("show",
  signature(object="SFFRead"),
  function(object){
    cat("Readname:\n")
    show(object@readname)
    cat("\n")
    cat("Bases:\n")
    show(object@bases)
    cat("\n")
    cat("Clip quality left:  ", object@clipQualityLeft, "\n", sep="")
    cat("Clip quality right: ", object@clipQualityRight, "\n", sep="")
    cat("Clip adapter left:  ", object@clipAdapterLeft, "\n", sep="")
    cat("Clip adapter right: ", object@clipAdapterRight, "\n", sep="")
    cat("\n")
    cat("Flowgram values: \n")
    lfv = length(object@flowgramValues)
    if(lfv > 20) {
      cat(object@flowgramValues[1:10], "...", object@flowgramValues[(lfv-10):lfv], "\n", sep=" ")
    } else {
      cat(object@flowgramValues)
    }
    cat("\n")
    cat("Flow indexes: \n")
    lfi = length(object@flowIndexes)
    if(lfi > 20) {
      cat(object@flowIndexes[1:10], "...", object@flowIndexes[(lfi-10):lfi], "\n", sep=" ")
    } else {
      cat(object@flowIndexes)
    }
    cat("\n")
    cat("Quality scores: \n")
    lqs = length(object@qualityScores)
    if(lqs > 20) {
      cat(object@qualityScores[1:10], "...", object@qualityScores[(lqs-10):lqs], "\n", sep=" ")
    } else {
      cat(object@qualityScores)
    }
  }
)
