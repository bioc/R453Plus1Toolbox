### Getter & Setter ###
setMethod("name", signature(x="SFFRead"),
  function(x) {
    return(x@name)
  }
)
setReplaceMethod("name", signature(x="SFFRead", value="character"),
  function(x, value) {
    x@name = value
    return(x)
  }
)

setMethod("read", signature(x="SFFRead"),
  function(x) {
    return(x@read)
  }
)
setReplaceMethod("read", signature(x="SFFRead", value="QualityScaledDNAStringSet"),
  function(x, value) {
    x@read = value
    return(x)
  }
)

setMethod("flowgramFormat", signature(x="SFFRead"),
  function(x) {
    return(x@flowgramFormat)
  }
)
setReplaceMethod("flowgramFormat", signature(x="SFFRead", value="numeric"),
  function(x, value) {
    x@flowgramFormat = value
    return(x)
  }
)

setMethod("keySequence", signature(x="SFFRead"),
  function(x) {
    return(x@keySequence)
  }
)
setReplaceMethod("keySequence", signature(x="SFFRead", value="character"),
  function(x, value) {
    x@keySequence = value
    return(x)
  }
)

setMethod("flowChars", signature(x="SFFRead"),
  function(x) {
    return(x@flowChars)
  }
)
setReplaceMethod("flowChars", signature(x="SFFRead", value="character"),
  function(x, value) {
    x@flowChars = value
    return(x)
  }
)

setMethod("clipQualityLeft", signature(x="SFFRead"),
  function(x) {
    return(x@clipQualityLeft)
  }
)
setReplaceMethod("clipQualityLeft", signature(x="SFFRead", value="numeric"),
  function(x, value) {
    x@clipQualityLeft = value
    return(x)
  }
)

setMethod("clipQualityRight", signature(x="SFFRead"),
  function(x) {
    return(x@clipQualityRight)
  }
)
setReplaceMethod("clipQualityRight", signature(x="SFFRead", value="numeric"),
  function(x, value) {
    x@clipQualityRight = value
    return(x)
  }
)

setMethod("clipAdapterLeft", signature(x="SFFRead"),
  function(x) {
    return(x@clipAdapterLeft)
  }
)
setReplaceMethod("clipAdapterLeft", signature(x="SFFRead", value="numeric"),
  function(x, value) {
    x@clipAdapterLeft = value
    return(x)
  }
)

setMethod("clipAdapterRight", signature(x="SFFRead"),
  function(x) {
    return(x@clipAdapterRight)
  }
)
setReplaceMethod("clipAdapterRight", signature(x="SFFRead", value="numeric"),
  function(x, value) {
    x@clipAdapterRight = value
    return(x)
  }
)

setMethod("flowgram", signature(x="SFFRead"),
  function(x) {
    return(x@flowgram)
  }
)
setReplaceMethod("flowgram", signature(x="SFFRead", value="numeric"),
  function(x, value) {
    x@flowgram = value
    return(object)
  }
)

setMethod("flowIndexes", signature(x="SFFRead"),
  function(x) {
    return(x@flowIndexes)
  }
)
setReplaceMethod("flowIndexes", signature(x="SFFRead", value="numeric"),
  function(x, value) {
    x@flowIndexes = value
    return(x)
  }
)

setMethod("show",
  signature(object="SFFRead"),
  function(object){
    cat("Readname:\n")
    show(object@name)
    cat("\n")
    cat("Flow chars:\n")
    show(object@flowChars)
    cat("\n")
    cat("Key sequence:\n")
    show(object@keySequence)
    cat("\n")
    cat("Read:\n")
    show(object@read)
    cat("\n")
    cat("Quality:\n")
    show(object@quality)
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
