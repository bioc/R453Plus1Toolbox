### Getter & Setter ###
setMethod("name", signature(object="SFFRead"),
  function(object) {
    return(object@name)
  }
)
setReplaceMethod("name", signature(object="SFFRead", value="character"),
  function(object, value) {
    object@name = value
    return(object)
  }
)

setMethod("read", signature(object="SFFRead"),
  function(object) {
    return(object@read)
  }
)
setReplaceMethod("read", signature(object="SFFRead", value="QualityScaledDNAStringSet"),
  function(object, value) {
    object@read = value
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

setMethod("flowgram", signature(object="SFFRead"),
  function(object) {
    return(object@flowgram)
  }
)
setReplaceMethod("flowgram", signature(object="SFFRead", value="numeric"),
  function(object, value) {
    object@flowgram = value
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

setMethod("show",
  signature(object="SFFRead"),
  function(object){
    cat("Readname:\n")
    show(object@name)
    cat("\n")
    cat("Read:\n")
    show(object@read)
    cat("\n")
    cat("Clip quality left:  ", object@clipQualityLeft, "\n", sep="")
    cat("Clip quality right: ", object@clipQualityRight, "\n", sep="")
    cat("Clip adapter left:  ", object@clipAdapterLeft, "\n", sep="")
    cat("Clip adapter right: ", object@clipAdapterRight, "\n", sep="")
    cat("\n")
    cat("Flowgram: \n")
    lfl = length(object@flowgram)
    if(lfl > 20) {
      cat(object@flowgram[1:10], "...", object@flowgram[(lfl-10):lfl], "\n", sep=" ")
    } else {
      cat(object@flowgram)
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
  }
)
