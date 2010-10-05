#################################################################################
########## METHODS TO CREATE A MAPPERSET ########################################
#################################################################################

setMethod("MapperSet",
    signature(dirs="character"),
    function(dirs, samplenames){

        ## catch some bad input
        ## first, check if all files in the given folders exist
        for(dir in dirs){
            file = file.path(dir, "mapping/454HCDiffs.txt")
	    if(!file.exists(file))
	    	stop("File not found: ", file)
	}

	## if no sample names are given, use the dirs
	if(!missing(samplenames)){
	    if(!(is.vector(samplenames, mode="character")))
	    	stop("invalid argument: please give a vector of characters as samplenames")
	    if(length(samplenames) != length(dirs))
		stop("invalid argument: number of sample names not consistent with the number of dirs")
    	}else
	    samplenames = dirs
        
        names(dirs) = samplenames

        message("Reading data ... ")
        mapperData1 = .readGSMAssayFeatureData(dirs, samplenames)
        mapperData2 = .readGSMPhenoData(dirs, samplenames)
        message("done")

        mapperSet = new("MapperSet",
            assay1=mapperData1[[1]],
	    assay2=mapperData1[[2]],
            assay3=mapperData1[[3]],
	    assay4=mapperData1[[4]],
            featureDf=mapperData1[[5]],
            featureMeta=mapperData1[[6]],
            phenoDf=mapperData2[[1]],
            phenoMeta=mapperData2[[2]],
            dirs=dirs
        )
        return(mapperSet)

    }
)

## extend the eSet initialize method

setMethod("initialize",
    signature(.Object="MapperSet"),
    function(.Object,
        assayData=assayDataNew("list", 
	    variantForwCount=assay1,
	    totalForwCount=assay2,
	    variantRevCount=assay3,
	    totalRevCount=assay4),
        featureData=new("AnnotatedDataFrame", data=featureDf,
            varMetadata=featureMeta),
        phenoData=new("AnnotatedDataFrame", data=phenoDf, varMetadata=phenoMeta),
        assay1=new("matrix"),
        assay2=new("matrix"),
        assay3=new("matrix"),
        assay4=new("matrix"),
        featureDf=new("data.frame"),
        featureMeta=new("data.frame"),
        phenoDf=new("data.frame"),
        phenoMeta=new("data.frame"),
        dirs="",
        ...)
    {
        .Object@dirs = dirs
        callNextMethod(.Object, assayData = assayData, featureData=featureData,
            phenoData=phenoData, ...)
    }
)

#################################################################################
########## IMPORT OF FEATURE AND PHENODATA (internal/private methods) ###########
#################################################################################

.readGSMAssayFeatureData <- function(dirs, samplenames) {

    firstDir = TRUE
    numSamples = length(dirs)
    sampleNAmes = NULL

    ## initialize assay and feature data (in case no file can be read)
    assayMutForward = matrix(NA, 1, numSamples)
    assayMutReverse = matrix(NA, 1, numSamples)
    assayForwardDepth = matrix(NA, 1, numSamples)
    assayReverseDepth = matrix(NA, 1, numSamples)
    featureData = data.frame(
        chr = NA,
        strand = NA,
        start = NA,
        end = NA,
        referenceBases = NA,
        variantBase = NA,
        regName = NA,
        knownSNP = NA,
        stringsAsFactors=FALSE
    )

    for(dir in dirs){
      

        message("... from ", dir)
	## read file
        file = file.path(dir, "mapping/454HCDiffs.txt")
       	text = readLines(file)

	## count reads with mutation
    	st = grep("Reads with Difference", text)
    	en = grep("Other Reads", text)
	
	##if the file is empty: only fill the assayData with zero values (no new featureData added)
	if(length(st) == 0 | length(en) == 0){
    	   assayMutForward = cbind(assayMutForward, matrix(0, nrow(assayMutForward),1))
    	   assayMutReverse = cbind(assayMutReverse, matrix(0, nrow(assayMutReverse),1))
    	   assayForwardDepth = cbind(assayForwardDepth, matrix(0, nrow(assayForwardDepth),1))
    	   assayReverseDepth = cbind(assayReverseDepth, matrix(0, nrow(assayReverseDepth),1))
	}
	else{

	    sampleNames = c(sampleNames, dir)

    	    countMutForward = vector(mode="numeric", length=length(st))
    	    countMutReverse = vector(mode="numeric", length=length(st))
    	    for(i in 1:length(st)) {
            	reads = text[(st[i]+2):(en[i]-1)]
	    	direction = grep("[A-Za-z0-9]", reads, value=TRUE)
	    	direction = gsub("[A-Za-z0-9 \\(\\)\\ ]", "", direction)
	    	direction = as.vector(sapply(direction, function(x) substring(x,1,1)))
            	countMutForward[i] = sum(direction == "+")
            	countMutReverse[i] = sum(direction == "-")
    	    }

	    ## count reads without mutation
    	    st = grep("Other Reads", text)
    	    en = grep(">chr", text)
    	    en = c(en[-1], length(text))
    	    countNoMutForward = vector(mode="numeric", length=length(st))
    	    countNoMutReverse = vector(mode="numeric", length=length(st))
    	    for(i in 1:length(st)) {
            	reads = text[(st[i]+1):(en[i]-1)]
	    	direction = grep("[A-Za-z0-9]", reads, value=TRUE)
	    	direction = gsub("[A-Za-z0-9 \\(\\)\\ ]", "", direction)
	    	direction = as.vector(sapply(direction, function(x) substring(x,1,1)))
            	countNoMutForward[i] = sum(direction == "+")
            	countNoMutReverse[i] = sum(direction == "-")
    	    }

	    ## read start and strand of reference sequence
    	    refLines = grep("^chr", text, value=TRUE)
    	    strands = vector(mode="character", length=length(refLines))
    	    refStarts = vector(mode="numeric", length=length(refLines))
    	    refSeqs = vector(mode="numeric", length=length(refLines))
    	    for(i in 1:length(refLines)) {
            	entries = refLines[i]
            	line = unlist(strsplit(entries, " +"))
            	strands[i] = substring(line[2], nchar(line[2]))
            	refStarts[i] = substring(line[2], 1, nchar(line[2])-1)
            	refSeqs[i] = gsub("-", "", line[3])
    	    }

	    ## construct eSet information containing assayData, featureData and phenoData	
    	    featLines = grep("^>chr", text, value=TRUE)
    	    con = textConnection(featLines)
    	    tab = read.csv(con, sep="\t", 
	    	col.names=c("Chr", "Start", "End", "RefNuc", "VarNuc", "TotalDepth", 
	    	    "VarFreq", "RefAA", "VarAA", "CodingFrame", "RegionName", "KnownSNPs"), 
	    	header=FALSE,
    	  	colClasses=c("character", "numeric", "numeric", "character", "character", 
	    	"numeric", "character", "character", "character", "character", "character", "character")
     	    )
    	    close(con)
	    chrs = gsub(">", "", tab[, "Chr"])
    	    tab[, "Chr"] = substr(chrs, 4, nchar(chrs))
    	    tab[,"mutReadsForward"] = countMutForward
    	    tab[,"mutReadsReverse"] = countMutReverse
    	    tab[,"forwardDepth"] = countMutForward + countNoMutForward
    	    tab[,"reverseDepth"] = countMutReverse + countNoMutReverse
    	    tab[,"Strand"] = strands

    	    ## case1 (first file): only initialize data
    	    if(firstDir){

    	    	## assay data
    	    	assayMutForward = matrix(countMutForward)
    	    	assayMutReverse = matrix(countMutReverse)
    	    	assayForwardDepth = matrix(countMutForward + countNoMutForward)
    	    	assayReverseDepth = matrix(countMutReverse + countNoMutReverse)

    	    	## feature data
    	        featureData = data.frame(
                   chr = tab[, "Chr"],
		   strand = tab[, "Strand"],
                   start = as.numeric(tab[, "Start"]),
                   end = as.numeric(tab[, "End"]),
                   referenceBases = tab[, "RefNuc"],
                   variantBase = tab[, "VarNuc"],
                   regName = tab[, "RegionName"],
                   knownSNP = tab[, "KnownSNPs"],
		   stringsAsFactors=FALSE
    	        )
	    	firstDir = FALSE
    	    }
    	    ## case2 (more than one file): merge new assay and feature data with older data
    	    else{

    	    	## assay data
    	    	assayMutForwardNew = matrix(countMutForward)
    	    	assayMutReverseNew = matrix(countMutReverse)
    	    	assayForwardDepthNew = matrix(countMutForward + countNoMutForward)
    	    	assayReverseDepthNew = matrix(countMutReverse + countNoMutReverse)

	    	## feature data
            	featureDataNew = data.frame(
                    chr = tab[, "Chr"],
		    strand = tab[, "Strand"],
                    start = tab[, "Start"],
                    end = tab[, "End"],
                    referenceBases = tab[, "RefNuc"],
                    variantBase = tab[, "VarNuc"],
                    regName = tab[, "RegionName"],
                    knownSNP = tab[, "KnownSNPs"],
		    stringsAsFactors=FALSE
	     	)
		## determine common rows between the old and the new feature data (then merge common data and add new entries)
		join <- function(x) do.call("paste", c(as.data.frame(x), sep = "\r"))
		idx1 = which(join(featureData) %in% join(featureDataNew))
		idx2 = which(join(featureDataNew) %in% join(featureData))
		idx1_not = setdiff(1:nrow(featureData), idx1)
		idx2_not = setdiff(1:nrow(featureDataNew), idx2)

	    	## extend assay data (by columns, initialization with 0)
    	    	assayMutForward = cbind(assayMutForward, matrix(0, nrow(assayMutForward),1))
    	    	assayMutReverse = cbind(assayMutReverse, matrix(0, nrow(assayMutReverse),1))
    	    	assayForwardDepth = cbind(assayForwardDepth, matrix(0, nrow(assayForwardDepth),1))
    	    	assayReverseDepth = cbind(assayReverseDepth, matrix(0, nrow(assayReverseDepth),1))

	    	## add entries for common mutations (assure, that both sets are in the same order (-> idx3))
		f1 = featureData[idx1, ]
	 	f2 = featureDataNew[idx2, ]
		a1 = assayMutForwardNew[idx2]
		a2 = assayMutReverseNew[idx2]
		a3 = assayForwardDepthNew[idx2]
		a4 = assayReverseDepthNew[idx2]
		idx3 = unlist(sapply(f1$start, function(x) which(f2$start == x)))

	        assayMutForward[idx1, ncol(assayMutForward)] = a1[idx3]
	        assayMutReverse[idx1, ncol(assayMutReverse)] = a2[idx3]
	        assayForwardDepth[idx1, ncol(assayForwardDepth)] = a3[idx3]
	        assayReverseDepth[idx1, ncol(assayReverseDepth)] = a4[idx3]

	        ## add entries for new mutations (by columns, insert 0 for older entries)
	   	## assay data slot 1
	    	aMutForwardNew = assayMutForwardNew[idx2_not]
	    	aMutForwardNew = cbind(matrix(0, length(aMutForwardNew), (ncol(assayMutForward)) - 1), aMutForwardNew)
	    	assayMutForward = rbind(assayMutForward, aMutForwardNew)
	    	## assay data slot 2
	    	aMutReverseNew = assayMutReverseNew[idx2_not]
	    	aMutReverseNew = cbind(matrix(0, length(aMutReverseNew), (ncol(assayMutReverse)) - 1), aMutReverseNew)
	    	assayMutReverse = rbind(assayMutReverse, aMutReverseNew)
	    	## assay data slot 3
	    	aForwardDepthNew = assayForwardDepthNew[idx2_not]
	    	aForwardDepthNew = cbind(matrix(0, length(aForwardDepthNew), (ncol(assayForwardDepth)) - 1), aForwardDepthNew)
	    	assayForwardDepth = rbind(assayForwardDepth, aForwardDepthNew)
	    	## assay data slot 4
	    	aReverseDepthNew = assayReverseDepthNew[idx2_not]
	    	aReverseDepthNew = cbind(matrix(0, length(aReverseDepthNew), (ncol(assayReverseDepth)) - 1), aReverseDepthNew)
	    	assayReverseDepth = rbind(assayReverseDepth, aReverseDepthNew)

	    	## extend feature data (by rows)
	    	fDataNew = featureDataNew[idx2_not, ]
	        featureData = rbind(featureData, fDataNew)
    	    }
        }
    }

    ## change style of varNuc: repeat "-" according to the number of bases involved in a deletion
    for(i in 1:nrow(featureData)){
	if(!is.na(featureData[i, "variantBase"]) & featureData[i, "variantBase"] == "-") 
	    featureData[i, "variantBase"] = paste(rep("-", nchar(featureData[i, "referenceBases"])), collapse="")
	if(!is.na(featureData[i, "referenceBases"]) & featureData[i, "referenceBases"] == "-") 
	    featureData[i, "referenceBases"] = paste(rep("-", nchar(featureData[i, "variantBase"])), collapse="")
    }

    ## distribute some new row and columnnames for feature and assay data
    rownames(featureData) = 1:nrow(featureData)
    rownames(assayMutForward) = rownames(featureData)
    rownames(assayMutReverse) = rownames(featureData)
    rownames(assayForwardDepth) = rownames(featureData)
    rownames(assayReverseDepth) = rownames(featureData)
    colnames(assayMutForward) = samplenames
    colnames(assayMutReverse) = samplenames
    colnames(assayForwardDepth) = samplenames
    colnames(assayReverseDepth) = samplenames

    ## feature meta
    featureMeta = data.frame(row.names=colnames(featureData),
        labelDescription=c("Chromosome","Strand", "Start", "End", "Reference nucleotide",
      	    "Variant nucleotide", "Region name", "Known SNPs"))

    return(
        list(
	    assayMutForward,
	    assayForwardDepth,
	    assayMutReverse,
	    assayReverseDepth,
	    featureData,
	    featureMeta
	)
    )


}

.readGSMPhenoData <- function(dirs, samplenames) {

    numSamples = length(dirs)

    ## pheno data
    accessionNumbers = rep(NA, numSamples)
    names(accessionNumbers) = dirs
    for(dir in dirs){

        file = file.path(dir, "mapping/454NewblerMetrics.txt")

        if(!file.exists(file))
            warning("File not found:", file, "; pheno data omitted")
        else{
            ## read file
            text = readLines(file)
            ## locate accession number
            an = grep(".sff", text, value=TRUE)
	    if(length(an) > 0){
            	an = an[1]
            	an = unlist(strsplit(an, split="/"))
            	an = an[length(an)]
            	an = unlist(strsplit(an, split=".sff"))
            	accessionNumbers[dir] = an[1]
	    }
	}
    }

    phenoData = matrix(accessionNumbers, numSamples, 1)
    phenoData = data.frame(accessionNumber=phenoData, row.names=samplenames)

    ## pheno meta
    phenoMeta = data.frame(row.names=colnames(phenoData),
        labelDescription=c("Accession number"))

    return(
        list(
            phenoData,
            phenoMeta
        )
    )
}




#################################################################################
########## GETTER AND SETTER METHODS ############################################
#################################################################################

## provide another subset method for samples and variants
## Todo: common method for AVASet & MapperSet
setMethod("subset",
    signature("MapperSet"),
          
    function(x, subset, dimension="variants"){

      object = x
      if(!missing(subset)){

            if(dimension=="variants"){
              ## catch some bad input
              if(is.vector(subset, mode="character")){
             	  if(!all(subset %in% rownames(fData(object)))){
                     stop("invalid argument: given variants not found")
             	  }
              }else{
                  if(is.vector(subset, mode="numeric")){
                      if(!all(subset %in% 1:nrow(fData(object)))){
                          stop("invalid argument: given indices out of range")
                      }
                  }else{
                      if(is.vector(subset, mode="logical")){
                          if(length(subset) != nrow(fData(object))){
                              stop("invalid argument: given subset is of wrong length")
                          }
                      }else
                          stop("invalid argument: given subset of incorrect type")
                  }
              }
              return(object[subset, ])
            }else{
              if(dimension=="samples"){
              ## catch some bad input
              if(is.vector(subset, mode="character")){
             	  if(!all(subset %in% rownames(pData(object)))){
                     stop("invalid argument: given samples not found")
             	  }
              }else{
                  if(is.vector(subset, mode="numeric")){
                      if(!all(subset %in% 1:nrow(pData(object)))){
                          stop("invalid argument: given indices out of range")
                      }
                  }else{
                      if(is.vector(subset, mode="logical")){
                          if(length(subset) != nrow(pData(object))){
                              stop("invalid argument: given subset is of wrong length")
                          }
                      }else
                          stop("invalid argument: given subset of incorrect type")
                  }
              }
                return(object[, subset])
              }else{
                stop("invalid argument: dimension can be either \"samples\" or \"variants\"")
              }
            }
          
        }else{
          warning("no subset given; MapperSet object left unchanged")
          return(object)
        }
    }
)

# this methods reads the file "454ReadStatus.txt" in the GSM project directory and returns its content in a dataframe
# the file contains information about the alignment of each read (chr, pos, strand, etc.)

setMethod("getReadStatus", signature(object="MapperSet"),
          function(object){
            dirs = dirs(object)
            readStatus = list()
            files = file.path(dirs, "mapping", "454ReadStatus.txt")
            for (i in 1:length(files)) {
              readStatus[[i]] = read.table(files[i], skip=1, header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
              names(readStatus)[i] = names(dirs)[i]
            }
            return(readStatus)
          })

# this version of assayData returns the filtered assayData of the MapperSet if the
# filter is set to a value > 0; otherwise it calls the original Biobase version

setMethod("assayData",
    signature(object="MapperSet"),
    function(object){
        if(sum(variantFilterPerc(object)) > 0){
            # filter the assayData
            filteredForwardCount=subset(object@assayData[[1]], rownames(object@assayData[[1]]) %in% variantFilter(object))
            filteredForwardDepth=subset(object@assayData[[2]], rownames(object@assayData[[2]]) %in% variantFilter(object))
            filteredReverseCount=subset(object@assayData[[3]], rownames(object@assayData[[3]]) %in% variantFilter(object))
            filteredReverseDepth=subset(object@assayData[[4]], rownames(object@assayData[[4]]) %in% variantFilter(object))
            return(assayDataNew("list",
                variantForwCount=filteredForwardCount,
                totalForwCount=filteredForwardDepth,
                variantRevCount=filteredReverseCount,
                totalRevCount=filteredReverseDepth))
        } else {
            callNextMethod(object)
        }
    }
)

# these versions of fetureData and fData return the filtered featureData of the MapperSet if the
# filter is set to a value > 0; otherwise it calls the original Biobase version

setMethod("fData",
    signature(object="MapperSet"),
    function(object){
        if(sum(variantFilterPerc(object)) > 0) {
            #filter the featureData
            return ( subset( object@featureData@data, rownames(object@featureData@data) %in% variantFilter(object) ) )
        } else {
            callNextMethod(object)
        }
    }
)

setMethod("featureData",
    signature(object="MapperSet"),
    function(object){
        if(sum(variantFilterPerc(object)) > 0) {
            #filter the featureData
            return ( subset( object@featureData, rownames(object@featureData@data) %in% variantFilter(object) ) )
        } else {
            callNextMethod(object)
        }
    }
)

setMethod("dirs",
    signature(object="MapperSet"),
    function(object){
        return(object@dirs)
    }
)

