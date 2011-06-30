setMethod("addRead", signature(object="SFFContainer", read="SFFRead"),
  function(object, read) {
    object@clipping[name(read), ] = c(clipQualityLeft(read), 
clipQualityRight(read), clipAdapterLeft(read), clipAdapterRight(read))
    object@flowgrams[[name(read)]] = flowgramValues(read)
    bases = bases(read)
    names(bases) = name(read)
    object@reads = c(object@reads, DNAStringSet(bases))
    object@flowIndexes[[name(read)]] = flowIndexes(read)
    return(object)
  }
)

setMethod("getRead", signature(object="SFFContainer", readname="character"),
  function(object, readname) {
    if(! readname %in% names(reads(object))) {
      stop(cat("There is no read with name ", readname, " in the SFFContainer.", sep=""))
    } else {
      cql = clipQualityLeft(object)[readname]
      cqr = clipQualityRight(object)[readname]
      cal = clipAdapterLeft(object)[readname]
      car = clipAdapterRight(object)[readname]
      names(cql) = names(cqr) = names(cal) = names(car) = NULL
      read = new("SFFRead", 
                 name=readname, 
                 bases=reads(object)[[readname]],
                 clipQualityLeft=cql,
                 clipQualityRight=cqr,
                 clipAdapterLeft=cal,
                 clipAdapterRight=car,
                 flowgramValues=flowgrams(object)[[readname]],
                 flowIndexes=flowIndexes(object)[[readname]]
      )
      return(read)
    }
  }
)

setMethod("name", signature(object="SFFContainer"),
  function(object) {
    return(object@name)
  }
)
setReplaceMethod("name", signature(object="SFFContainer", value="character"),
  function(object, value) {
    object@name = value
    return(object)
  }
)

setMethod("flowgramFormat", signature(object="SFFContainer"),
  function(object) {
    return(object@flowgramFormat)
  }
)
setReplaceMethod("flowgramFormat", signature(object="SFFContainer", value="numeric"),
  function(object, value) {
    object@flowgramFormat = value
    return(object)
  }
)

setMethod("keySequence", signature(object="SFFContainer"),
  function(object) {
    return(object@keySequence)
  }
)
setReplaceMethod("keySequence", signature(object="SFFContainer", value="character"),
  function(object, value) {
    object@keySequence = value
    return(object)
  }
)

setMethod("flowChars", signature(object="SFFContainer"),
  function(object) {
    return(object@flowChars)
  }
)
setReplaceMethod("flowChars", signature(object="SFFContainer", value="character"),
  function(object, value) {
    object@flowChars = value
    return(object)
  }
)

setMethod("clipQualityLeft", signature(object="SFFContainer"),
  function(object) {
    return(object@clipQualityLeft)
  }
)
setReplaceMethod("clipQualityLeft", signature(object="SFFContainer", value="numeric"),
  function(object, value) {
    object@clipQualityLeft = value
    return(object)
  }
)

setMethod("clipQualityRight", signature(object="SFFContainer"),
  function(object) {
    return(object@clipQualityRight)
  }
)
setReplaceMethod("clipQualityRight", signature(object="SFFContainer", value="numeric"),
  function(object, value) {
    object@clipQualityRight = value
    return(object)
  }
)

setMethod("clipAdapterLeft", signature(object="SFFContainer"),
  function(object) {
    return(object@clipAdapterLeft)
  }
)
setReplaceMethod("clipAdapterLeft", signature(object="SFFContainer", value="numeric"),
  function(object, value) {
    object@clipAdapterLeft = value
    return(object)
  }
)

setMethod("clipAdapterRight", signature(object="SFFContainer"),
  function(object) {
    return(object@clipAdapterRight)
  }
)
setReplaceMethod("clipAdapterRight", signature(object="SFFContainer", value="numeric"),
  function(object, value) {
    object@clipAdapterRight = value
    return(object)
  }
)

setMethod("flowgrams", signature(object="SFFContainer"),
  function(object) {
    return(object@flowgrams)
  }
)
setReplaceMethod("flowgrams", signature(object="SFFContainer", value="list"),
  function(object, value) {
    object@flowgrams = value
    return(object)
  }
)

setMethod("flowIndexes", signature(object="SFFContainer"),
  function(object) {
    return(object@flowIndexes)
  }
)
setReplaceMethod("flowIndexes", signature(object="SFFContainer", value="list"),
  function(object, value) {
    object@flowIndexes = value
    return(object)
  }
)

setMethod("reads", signature(object="SFFContainer"),
  function(object) {
    return(object@reads)
  }
)
setReplaceMethod("reads", signature(object="SFFContainer", value="QualityScaledDNAStringSet"),
  function(object, value) {
    object@reads = value
    return(object)
  }
)

setMethod("show",
  signature(object="SFFContainer"),
  function(object){
    cat("Name: \n")
    cat(object@name, "\n")
    cat("\n")
    cat("Reads: \n")
    show(object@reads)
    cat("\n")
  }
)

setMethod("[",
  signature("SFFContainer"),
  function(x, i, ...){
    if(is.numeric(i)) {
      sffNew = new("SFFContainer", 
                   name=name(x),
                   flowgramFormat=flowgramFormat(x),
                   flowChars=flowChars(x),
                   keySequence=keySequence(x),
                   clipQualityLeft=clipQualityLeft(x)[i],
                   clipQualityRight=clipQualityRight(x)[i],
                   clipAdapterLeft=clipAdapterLeft(x)[i],
                   clipAdapterRight=clipAdapterRight(x)[i],
                   flowgrams=flowgrams(x)[i],
                   flowIndexes=flowIndexes(x)[i],
                   reads=reads(x)[i]
                  )
    } else {
      reads = reads(x)[names(reads) %in% i]
      sffNew = new("SFFContainer", 
                   name=name(x),
                   flowgramFormat=flowgramFormat(x),
                   flowChars=flowChars(x),
                   keySequence=keySequence(x),
                   clipQualityLeft=clipQualityLeft(x)[i],
                   clipQualityRight=clipQualityRight(x)[i],
                   clipAdapterLeft=clipAdapterLeft(x)[i],
                   clipAdapterRight=clipAdapterRight(x)[i],
                   flowgrams=flowgrams(x)[i],
                   flowIndexes=flowIndexes(x)[i],
                   reads=reads
                  )
    }
    return(sffNew)
  }
)

setMethod("sff2fastq", 
  signature(object="SFFContainer"), 
  function(object, outdir, fname){
    if(missing(outdir)) {
      outdir = getwd()
    }
    outdir = gsub(paste(.Platform$file.sep, "+$", sep=""), "", outdir)
    if(missing(fname)) {
      fname = name(object)
      fname = gsub("(\\.sff)?$", ".fastq", fname)
    }
    fp = file.path(outdir, fname)
    write.XStringSet(reads(object), filepath=fp, format="fastq", qualities=quality(reads(object)))
    cat("Written file", fp, "\n")
  }
)
