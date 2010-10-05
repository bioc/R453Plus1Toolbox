##################################################################################################
## getVariantPercentages returns the variant coverage (forward and/or reverse) in percent
## optionally: filters variants whose coverage is higher or equal than a given threshold
## output: data frame with the (filtered) coverage for each variant
##################################################################################################

.getVariantPercentages <- function(object, direction="both"){

    numFilters = length(filter)

    if(!is.character(direction) || !(direction %in% c("forward", "reverse", "both")))
        stop("invalid argument: please specify \"forward\", \"reverse\" or \"both\" as direction")

    if(direction == "both")
    	percs = (assayData(object)[[1]] + assayData(object)[[3]]) /
             	(assayData(object)[[2]] + assayData(object)[[4]])
    if(direction == "forward")
    	percs = assayData(object)[[1]] / assayData(object)[[2]]
    if(direction == "reverse")
    	percs = assayData(object)[[3]] / assayData(object)[[4]]

#    percs[is.na(percs)] = 0
    return(as.matrix(percs))
}


######################################################################################################
## The method setVariantFilter sets the filter for variant coverage to a given value between 0 and 1
## additionally it sets the data frame filter for variants for future access to feature and assayData.
## It returns an updated AVAset / MapperSet
######################################################################################################

.setVariantFilter <- function(object, filter = 0){

    numFilters = length(filter)

    ## catch bad filter input (only filter values between 0 and 1 allowed)
    if(!is.numeric(filter) || !(numFilters %in% c(1,2)))
        stop("invalid argument: please specify one or two numeric filter values between 0 and 1")
    if(any(filter < 0) || any(filter > 1))
        stop("argument out of range: please specify one or two numeric filter values between 0 and 1")

    ## reset filter
    variantFilterPerc(object) = 0

    numVars = nrow(fData(object))
    numSamples = nrow(pData(object))
    variantsPercFiltered = vector("list", 2)
    for(i in 1:numFilters){

	f = filter[i]
	## if one filter value is 0: take all variants without filtering
        if(f == 0){
	    variantsPercFiltered[[i]] = matrix(1, numVars, numSamples)
	    rownames(variantsPercFiltered[[i]]) = rownames(assayData(object)[[1]])
	    colnames(variantsPercFiltered[[i]]) = colnames(assayData(object)[[1]])
        } else {

            ## build the filter:
	    ## if just one filter value specified: apply filter to the sum of both directions (forward and reverse)
	    if(numFilters == 1)
            	variantsPerc = (assayData(object)[[1]] + assayData(object)[[3]]) /
                               (assayData(object)[[2]] + assayData(object)[[4]])
    	    ## if two filter values specified: apply filter to each direction (forward / reverse)
	    else
		if(i == 1)
                    variantsPerc = assayData(object)[[1]] / assayData(object)[[2]]
		else
                    variantsPerc = assayData(object)[[3]] / assayData(object)[[4]]

	    ## remove NA-values and mark those variants, whose percentage is larger than the filter (perc >= filter in forw.&rev. for at least one sample)
            variantsPerc[is.na(variantsPerc)] = 0
            variantsPerc[variantsPerc < f] = 0
            variantsPerc[variantsPerc >= f] = 1
	    variantsPercFiltered[[i]] = variantsPerc
        }
    }
    ## if just one filter value specified: just return the (first) filter
    if(numFilters == 1){
        names(filter) = c("forward & reverse")
        variantFilterPerc(object) = filter
	validVariants = apply(variantsPercFiltered[[1]] == 1, 1, any)
    	variantFilter(object) = names(validVariants[validVariants])
        message("combined forward & reverse filter set to ", filter)
    ## if two filter values specified: return the variant names that meet both filter requirements
    }else{
        names(filter) = c("forward", "reverse")
        variantFilterPerc(object) = filter
   	## the filter finally consists of the variant names that meet the filter requirements (perc >= filter in forw.&rev. for at least one sample)
	validVariants = apply(variantsPercFiltered[[1]] & variantsPercFiltered[[2]], 1, any)
    	variantFilter(object) = names(validVariants[validVariants])
        message("forward filter set to ", filter[1])
        message("reverse filter set to ", filter[2])
    }

    return(object)
}


setMethod("getVariantPercentages",
    signature(object="AVASet"),
    .getVariantPercentages)

setMethod("getVariantPercentages",
    signature(object="MapperSet"),
    .getVariantPercentages)

setMethod("setVariantFilter",
    signature(object="AVASet"),
	.setVariantFilter)	

setMethod("setVariantFilter",
    signature(object="MapperSet"),
	.setVariantFilter)


## getter and setter for the AVASet

setMethod("variantFilterPerc",
    signature(object="AVASet"),
    function(object){
        return(object@variantFilterPerc)
    }
)

setMethod("variantFilter",
    signature(object="AVASet"),
    function(object){
        return(object@variantFilter)
    }
)
setReplaceMethod("variantFilterPerc",
    signature(object="AVASet", value="numeric"),
    function(object, value){
        object@variantFilterPerc=value
        return(object)
    }
)

setReplaceMethod("variantFilter",
    signature(object="AVASet", value = "character"),
    function(object, value){
        object@variantFilter=value
        return(object)
    }
)

## getter and setter for the MapperSet

setMethod("variantFilterPerc",
    signature(object="MapperSet"),
    function(object){
        return(object@variantFilterPerc)
    }
)

setMethod("variantFilter",
    signature(object="MapperSet"),
    function(object){
        return(object@variantFilter)
    }
)

setReplaceMethod("variantFilterPerc",
    signature(object="MapperSet", value="numeric"),
    function(object, value){
        object@variantFilterPerc=value
        return(object)
    }
)

setReplaceMethod("variantFilter",
    signature(object="MapperSet", value = "character"),
    function(object, value){
        object@variantFilter=value
        return(object)
    }
)

