setMethod("addRead", signature(x="SFFContainer", read="SFFRead"),
  function(x, read) {
    if(name(read) %in% names(reads(x))) {
      stop("A read with name ", name(read), " already exists.")
    } else if (flowgramFormat(read) != flowgramFormat(x)) {
      stop("Flowgram formats did not match: ", flowgramFormat(read), " vs. ", flowgramFormat(x))
    } else if (flowChars(read) != flowChars(x)) {
      stop("Flow chars did not match: ", flowChars(read), " vs. ", flowChars(x))
    } else if (keySequence(read) != keySequence(x)) {
      stop("Key sequences did not match: ", keySequence(read), " vs. ", keySequence(x))
    } else {
      x@clipQualityLeft[name(read)] = clipQualityLeft(read)
      x@clipQualityRight[name(read)] = clipQualityRight(read)
      x@clipAdapterLeft[name(read)] = clipAdapterLeft(read)
      x@clipAdapterRight[name(read)] = clipAdapterRight(read)
      x@flowgrams[[name(read)]] = flowgram(read)
      x@flowIndexes[[name(read)]] = flowIndexes(read)
      
      pq = PhredQuality(quality(read))
      names(pq) = name(read)
      newRead = QualityScaledDNAStringSet(read(read), pq)
      names(newRead) = name(read)
      
      origReads = reads(x)
      
      # didn't work
      # newobj = append(origReads, newRead)
      # x@reads = newobj
      # workaround
      s1 = append(as(newRead, "XStringSet"), as(origReads, "XStringSet"))
      s2 = as(append(newRead@quality, origReads@quality), class(origReads@quality))
      x@reads = QualityScaledDNAStringSet(s1, s2)
      return(x)
    }
  }
)

setMethod("getRead", signature(x="SFFContainer", readname="character"),
  function(x, readname) {
    if(! readname %in% names(reads(x))) {
      stop(cat("There is no read with name ", readname, " in the SFFContainer.", sep=""))
    } else {
      cql = clipQualityLeft(x)[readname]
      cqr = clipQualityRight(x)[readname]
      cal = clipAdapterLeft(x)[readname]
      car = clipAdapterRight(x)[readname]
      names(cql) = names(cqr) = names(cal) = names(car) = NULL
      read = new("SFFRead", 
                 name=readname,
                 flowgramFormat=flowgramFormat(x),
                 flowChars=flowChars(x),
                 keySequence=keySequence(x),
                 clipQualityLeft=cql,
                 clipQualityRight=cqr,
                 clipAdapterLeft=cal,
                 clipAdapterRight=car,
                 flowgram=flowgrams(x)[[readname]],
                 flowIndexes=flowIndexes(x)[[readname]],
                 read=reads(x)[[readname]],
                 quality=quality(reads(x))[[readname]]
      )
      return(read)
    }
  }
)

setMethod("name", signature(x="SFFContainer"),
  function(x) {
    return(x@name)
  }
)
setReplaceMethod("name", signature(x="SFFContainer", value="character"),
  function(x, value) {
    x@name = value
    return(x)
  }
)

setMethod("flowgramFormat", signature(x="SFFContainer"),
  function(x) {
    return(x@flowgramFormat)
  }
)
setReplaceMethod("flowgramFormat", signature(x="SFFContainer", value="numeric"),
  function(x, value) {
    x@flowgramFormat = value
    return(x)
  }
)

setMethod("keySequence", signature(x="SFFContainer"),
  function(x) {
    return(x@keySequence)
  }
)
setReplaceMethod("keySequence", signature(x="SFFContainer", value="character"),
  function(x, value) {
    x@keySequence = value
    return(x)
  }
)

setMethod("flowChars", signature(x="SFFContainer"),
  function(x) {
    return(x@flowChars)
  }
)
setReplaceMethod("flowChars", signature(x="SFFContainer", value="character"),
  function(x, value) {
    x@flowChars = value
    return(x)
  }
)

setMethod("clipQualityLeft", signature(x="SFFContainer"),
  function(x) {
    return(x@clipQualityLeft)
  }
)
setReplaceMethod("clipQualityLeft", signature(x="SFFContainer", value="numeric"),
  function(x, value) {
    x@clipQualityLeft = value
    return(x)
  }
)

setMethod("clipQualityRight", signature(x="SFFContainer"),
  function(x) {
    return(x@clipQualityRight)
  }
)
setReplaceMethod("clipQualityRight", signature(x="SFFContainer", value="numeric"),
  function(x, value) {
    x@clipQualityRight = value
    return(x)
  }
)

setMethod("clipAdapterLeft", signature(x="SFFContainer"),
  function(x) {
    return(x@clipAdapterLeft)
  }
)
setReplaceMethod("clipAdapterLeft", signature(x="SFFContainer", value="numeric"),
  function(x, value) {
    x@clipAdapterLeft = value
    return(x)
  }
)

setMethod("clipAdapterRight", signature(x="SFFContainer"),
  function(x) {
    return(x@clipAdapterRight)
  }
)
setReplaceMethod("clipAdapterRight", signature(x="SFFContainer", value="numeric"),
  function(x, value) {
    x@clipAdapterRight = value
    return(x)
  }
)

setMethod("flowgrams", signature(x="SFFContainer"),
  function(x) {
    return(x@flowgrams)
  }
)
setReplaceMethod("flowgrams", signature(x="SFFContainer", value="list"),
  function(x, value) {
    x@flowgrams = value
    return(x)
  }
)

setMethod("flowIndexes", signature(x="SFFContainer"),
  function(x) {
    return(x@flowIndexes)
  }
)
setReplaceMethod("flowIndexes", signature(x="SFFContainer", value="list"),
  function(x, value) {
    x@flowIndexes = value
    return(x)
  }
)

setMethod("reads", signature(x="SFFContainer"),
  function(x) {
    return(x@reads)
  }
)
setReplaceMethod("reads", signature(x="SFFContainer", value="QualityScaledDNAStringSet"),
  function(x, value) {
    x@reads = value
    return(x)
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
  signature(x="SFFContainer"), 
  function(x, outdir, fname){
    if(missing(outdir)) {
      outdir = getwd()
    }
    outdir = gsub(paste(.Platform$file.sep, "+$", sep=""), "", outdir)
    if(missing(fname)) {
      fname = name(x)
      fname = gsub("(\\.sff)?$", ".fastq", fname)
    }
    fp = file.path(outdir, fname)
    writeXStringSet(reads(x), filepath=fp, format="fastq", qualities=quality(reads(x)))
    cat("Written file", fp, "\n")
  }
)
