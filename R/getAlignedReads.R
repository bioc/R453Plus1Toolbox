.getAlignedReads = function(object, amplicons, dir){

    ## catch some bad input
    if(!missing(amplicons)){
    	if(!is.vector(amplicons, mode="character"))
	    stop("invalid argument: please specify your amplicon names as a charcter vector (see fDataAmp() for details)")
	if(!all(amplicons %in% rownames(fDataAmp(object))))
	    stop("invalid argument: (some of) your amplicon names were not found in the AVA set (see fDataAmp() for details)")
    }else
	amplicons = rownames(fDataAmp(object))

    ## directory of alignment file is either given directly by a function argument or is read from the AVASet object
    if(!missing(dir)){
    	if(!is.character(dir))
	    stop("invalid argument: directory has to be of type character")
	rootDir = file.path(dir, "Amplicons/Results/Align")
    }else
	rootDir = dirs(object)["align"]
    if(!file.exists(rootDir))
        stop("directory not found: ", rootDir)


    sampleIDs = pData(object)$SampleID
    refSeqs = fDataAmp(object)[amplicons, "referenceSeqID"]
    allReads = vector(mode="character")

    ## read individual file for each sample and every amplicon
    cat(length(amplicons), " amplicons were found\n")
    cat("Reading aligned reads from ", rootDir, " ... ")
    for(sample in sampleIDs){

    	for(ref in refSeqs){
	
	    file = paste(sample, "_vs_", ref, ".readAlign.txt", sep="")
	    file = file.path(rootDir, sample, ref, file)

	    if(file.exists(file)){

              ## read aligned reads; format: two rows per read, one name and one sequence
              df = read.table(file, fill=TRUE, sep="\t", skip=2, stringsAsFactors=FALSE, col.names=c("accessionNumber", "read", "tmp"))
              accessionNumbers = df$accessionNumber[seq(1, nrow(df), 2)] 		## remove entries from rows with sequences
              reads = df$read[df$read != ""]                              	## remove (empty) entries from rows with read names

              ## extract accession numbers from other information in that row
              accessionNumbers = strsplit(accessionNumbers, split = " /")
              accessionNumbers = sapply(accessionNumbers, function(x) substr(x[1], 2, nchar(x[1])))	    
              
              ## remove filling characters ("-") from read sequence (for example: reads have the form "TCTTCTTCAC-AGGTGCTTT-CAA-G-AAC-AGG--G")
              reads = gsub("-","",reads)

              names(reads) = accessionNumbers
              allReads = c(allReads, reads)
            }
            else
              warning("file not found: ", file)
	}
    }
    cat("done\n")
    cat(length(allReads), " aligned reads were found\n")
    
    ## return the reads as a DNAStringSet
    return(DNAStringSet(allReads, use.names=TRUE))

}

setMethod("getAlignedReads",
    signature(object="AVASet"),
    .getAlignedReads)
