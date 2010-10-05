## length, names are the same for seqsC1, commonBpsC1, commonAlignC1, alignedReadsC1 (and C2 respectively)

setMethod("length", signature(x="Breakpoints"),
    function(x) {
        length(commonBpsC1(x))
    }
)

setMethod("names", signature(x="Breakpoints"),
    function(x) {
        names(commonBpsC1(x))
    }
)

setReplaceMethod("names", signature(x="Breakpoints", value="ANY"),
    function(x, value) {
        names(seqsC1(x)) = value
        names(seqsC2(x)) = value
        names(commonBpsC1(x)) = value
        names(commonBpsC2(x)) = value
        names(commonAlignC1(x)) = value
        names(commonAlignC2(x)) = value
        names(alignedReadsC1(x)) = value
        names(alignedReadsC2(x)) = value
        return(x)
    }
)

setMethod("summary", signature(object="Breakpoints"),
    function(object) {

	if(length(object) > 0){
          chrA=sapply(commonBpsC1(object), function(x) as.character(x[1, "chr"]))
          chrB=sapply(commonBpsC1(object), function(x) as.character(x[2, "chr"]))
          bpAC1 = sapply(commonBpsC1(object), function(x) x[1, "breakpoint"])
          bpBC1 = sapply(commonBpsC1(object), function(x) x[2, "breakpoint"])
          bpAC2 = sapply(commonBpsC2(object), function(x) if(nrow(x) > 0) x[1, "breakpoint"] else return(NA))
          bpBC2 = sapply(commonBpsC2(object), function(x) if(nrow(x) > 0) x[2, "breakpoint"] else return(NA))

          # get chromosomes from Case 2 temporarly
          chrAC2=sapply(commonBpsC2(object), function(x) if(nrow(x) > 0) as.character(x[1, "chr"]) else return(NA))
          chrBC2=sapply(commonBpsC2(object), function(x) if(nrow(x) > 0) as.character(x[2, "chr"]) else return(NA))

          ind = !is.na(chrAC2)
          # chrA from case 2 is not chrA from case 1
          ind2 = chrA[ind] != chrAC2[ind]
          bpTmp = bpAC2[ind][ind2]
          bpAC2[ind][ind2] = bpBC2[ind][ind2]
          bpBC2[ind][ind2] = bpTmp
          # inversion
          ind = which(chrA == chrB & !is.na(chrAC2))
          for (i in ind) {
            if (abs(bpAC1[i] - bpAC2[i]) > abs(bpAC1[i] - bpBC2[i])) {
              bpTmp = bpAC2[i]
              bpAC2[i] = bpBC2[i]
              bpBC2[i] = bpTmp
            }
          }
          rm(chrAC2, chrBC2)


          # Get number of reads.
          firstCase = sapply(alignedReadsC1(object), function(x) {
            t = table(strand(x))
            if (length(t) == 1) {
              return(paste(t, "0", sep="/"))
            } else if (length(t) == 2) {
              return(paste(t[1], t[2], sep="/"))
            } else {
              stop("More than two strands!")
            }
          })
          secondCase = sapply(alignedReadsC2(object), function(x) {
            if (length(x) > 0) {
              t = table(strand(x))
              if (length(t) == 1) {
                return(paste(t, "0", sep="/"))
              } else if (length(t) == 2) {
                return(paste(t[1], t[2], sep="/"))
              } else {
                stop("More than two strands!")
              }
            } else {
              return("0/0")
            }
          })
          total = sapply(seqsC1(object), function(x) nrow(x)) + sapply(seqsC2(object), function(x) nrow(x))
          
          df = data.frame(
            ChrA=chrA,
            ChrB=chrB,
            BpACase1=bpAC1,
            BpBCase1=bpBC1,
            BpACase2=bpAC2,
            BpBCase2=bpBC2,
            NoReadsCase1=firstCase,
            NoReadsCase2=secondCase,
            NoReadsTotal=total)
        
          return(df)
        }else{
          warning("The given object is empty")
          return(data.frame())
        }
    }
)


setMethod(show, signature(object="Breakpoints"),
    function(object) {
        n = min(6, length(object))
 	if(n > 0){
            df = data.frame(
            	Size=sapply(seqsC1(object), function(x) nrow(x)) + sapply(seqsC2(object), function(x) nrow(x)),
            	ChrA=sapply(commonBpsC1(object), function(x) as.character(x[1, "chr"])),
            	ChrB=sapply(commonBpsC1(object), function(x) as.character(x[2, "chr"])),
            	stringsAsFactors=FALSE
            )
            if (n == 1) {
            	message(paste("Collection of 1 breakpoint:"))
            } else {
            	message(paste("Collection of", length(object), "breakpoints:"))
            }
            print(df[1:n,])
            if (length(object) > n) {
            	message("...")
            }
	}else
	    message("No consensus breakpoints detected: Object is empty")
    }
)

setMethod(table, signature(...="Breakpoints"),
    function(...) {
      if(length(...) > 0){
        size = sapply(seqsC1(...), function(x) nrow(x)) + sapply(seqsC2(...), function(x) nrow(x))
        return(table(size))
      }else{
        stop("The given object is empty")
      }
    }
)


#########
## GETTER
#########

setMethod("[",
    signature("Breakpoints"),
    function(x, i, ...){
        if((is.numeric(i) && all(is.element(i, -length(x):length(x))))  || is.element(i, names(x))) {
	    seqsC1(x) = seqsC1(x)[i]
	    seqsC2(x) = seqsC2(x)[i]
	    commonBpsC1(x) = commonBpsC1(x)[i]
	    commonBpsC2(x) = commonBpsC2(x)[i]
	    commonAlignC1(x) = commonAlignC1(x)[i]
	    commonAlignC2(x) = commonAlignC2(x)[i]
	    alignedReadsC1(x) = alignedReadsC1(x)[i]
	    alignedReadsC2(x) = alignedReadsC2(x)[i]
	    return(x)
	}else
            stop("Index invalid or out of range")
    }
)

setMethod("seqsC1", signature(object="Breakpoints"),
    function(object) {
        return(object@seqsC1)
    }
)
setMethod("seqsC2", signature(object="Breakpoints"),
    function(object) {
        return(object@seqsC2)
    }
)

setMethod("commonBpsC1", signature(object="Breakpoints"),
    function(object) {
        return(object@commonBpsC1)
    }
)
setMethod("commonBpsC2", signature(object="Breakpoints"),
    function(object) {
        return(object@commonBpsC2)
    }
)

setMethod("commonAlignC1", signature(object="Breakpoints"),
    function(object) {
        return(object@commonAlignC1)
    }
)
setMethod("commonAlignC2", signature(object="Breakpoints"),
    function(object) {
        return(object@commonAlignC2)
    }
)

setMethod("alignedReadsC1", signature(object="Breakpoints"),
    function(object) {
        return(object@alignedReadsC1)
    }
)
setMethod("alignedReadsC2", signature(object="Breakpoints"),
    function(object) {
        return(object@alignedReadsC2)
    }
)

#########
## SETTER
#########

setReplaceMethod("seqsC1",
    signature=signature(object="Breakpoints", value="list"),
    function(object, value) {
        object@seqsC1 = value
        return(object)
    }
)
setReplaceMethod("seqsC2",
    signature=signature(object="Breakpoints", value="list"),
    function(object, value) {
        object@seqsC2 = value
        return(object)
    }
)

setReplaceMethod("commonBpsC1",
    signature=signature(object="Breakpoints", value="list"),
    function(object, value) {
        object@commonBpsC1 = value
        return(object)
    }
)
setReplaceMethod("commonBpsC2",
    signature=signature(object="Breakpoints", value="list"),
    function(object, value) {
        object@commonBpsC2 = value
        return(object)
    }
)

setReplaceMethod("commonAlignC1",
    signature=signature(object="Breakpoints", value="list"),
    function(object, value) {
        object@commonAlignC1 = value
        return(object)
    }
)
setReplaceMethod("commonAlignC2",
    signature=signature(object="Breakpoints", value="list"),
    function(object, value) {
        object@commonAlignC2 = value
        return(object)
    }
)

setReplaceMethod("alignedReadsC1",
    signature=signature(object="Breakpoints", value="list"),
    function(object, value) {
        object@alignedReadsC1 = value
        return(object)
    }
)
setReplaceMethod("alignedReadsC2",
    signature=signature(object="Breakpoints", value="list"),
    function(object, value) {
        object@alignedReadsC2 = value
        return(object)
    }
)
