#################################################################################
########## METHODS TO CREATE AN AVASET ##########################################
#################################################################################

# AVASet reads all relevant information (sample-, variant- and amplicon
# data)
# from the (sub-) directories of the project
# output: AVASet, object from a subclass of the Biobase eSet and extended by
# two slots for amplicon data (assay- and feature data) and reference sequences

setMethod("AVASet",
    signature(dirname="character"),

    function(dirname){

        # check if the given dirname leads to all relevant files and directories
	dir_root = file.path(dirname, "Amplicons")
	dir_results = file.path(dir_root,"Results")
	dir_projectDef = file.path(dir_root,"ProjectDef")
	dir_variants = file.path(dir_results, "Variants")
	dir_align = file.path(dir_results, "Align")
	dirs = c(dir_results, dir_projectDef, dir_variants, dir_align)
	names(dirs) = c("results", "projectDef", "variants", "align")

        if(!file.exists(dirname)
            | !file.exists(dir_root) 
            | !file.exists(dir_results) 
            | !file.exists(dir_projectDef)
            | !file.exists(dir_variants)
            | !file.exists(dir_align)
            | !file.exists(file.path(dir_variants, "currentVariantDefs.txt"))
            | !file.exists(file.path(dir_projectDef, "ampliconsProject.txt"))
        )
            stop("No AVA project found under given directory or files in the directories are missing")

	# if all files were found, read sample, variant and amplicon data
        message("Reading sample data ... ", appendLF=FALSE)
        sampleData = readSampleData(dir_projectDef)

        message("done\nReading reference sequences ... ", appendLF=FALSE)
        refSeqData = readReferenceSequences(dir_projectDef)

        #reading all variant data requires the sample data and the reference
        # sequences
        message("done\nReading variant data ... ", appendLF=FALSE)
        varData = readVariants(dir_variants, sampleData[[1]], sread(refSeqData))

        #reading all amplicon data requires the sample data
        message("done\nReading amplicon data ... ", appendLF=FALSE)
        ampData = readAmplicons(dir_projectDef, dir_align, sampleData[[1]])
        message("done")

        avaSet = new("AVASet",
            variantForwCount=varData[[1]], 
            totalForwCount=varData[[2]], 
            variantRevCount=varData[[3]],
	    totalRevCount=varData[[4]], 
            forwCount=ampData[[1]],
            revCount=ampData[[2]],
            featureDf=varData[[5]], 
            featureMeta=varData[[6]], 
            phenoDf=sampleData[[1]], 
            phenoMeta=sampleData[[2]],
            featureDfAmp=ampData[[3]], 
            featureMetaAmp=ampData[[4]],
            refSeqs=refSeqData,
	    dirs = dirs
        )
        return(avaSet)
 
    }
)


# extend the eSet initialize method by the two new amplicon slots

setMethod("initialize",
    signature(.Object="AVASet"),
    function(.Object,
        assayData=assayDataNew("list", variantForwCount=variantForwCount,
            totalForwCount=totalForwCount,
            variantRevCount=variantRevCount,
            totalRevCount=totalRevCount),
        assayDataAmp=assayDataNew("list",
            forwCount=forwCount,
            revCount=revCount),
        featureData=new("AnnotatedDataFrame", data=featureDf,
            varMetadata=featureMeta),
        phenoData=new("AnnotatedDataFrame", data=phenoDf, varMetadata=phenoMeta),
        featureDataAmp=new("AnnotatedDataFrame", data=featureDfAmp,
            varMetadata=featureMetaAmp),
        variantForwCount=new("matrix"),
        variantRevCount=new("matrix"),
        totalForwCount=new("matrix"),
        totalRevCount=new("matrix"),
        forwCount=new("matrix"),
        revCount=new("matrix"),
        featureDf=new("data.frame"),
        featureMeta=new("data.frame"),
        phenoDf=new("data.frame"),
        phenoMeta=new("data.frame"),
        featureDfAmp=new("data.frame"),
        featureMetaAmp=new("data.frame"),
        refSeqs=new("AlignedRead"),
	dirs="",
        ...)
    {
        .Object@assayDataAmp = assayDataAmp
        .Object@featureDataAmp= featureDataAmp
        .Object@referenceSequences = refSeqs
	.Object@dirs = dirs
        callNextMethod(.Object, assayData = assayData, featureData=featureData,
            phenoData=phenoData, ...)
    }
)


#################################################################################
########## IMPORT OF FEATURE AND PHENODATA (internal methods) ###################
#################################################################################


# readSampleData runs through the subdirectory ".../Amplicons/Results/Variants"
# and reads the sample names and group infos

setMethod("readSampleData",
    signature(dir_projectDef="character"),
    function(dir_projectDef){

    # read the sample data (internal id):
    # file has no data frame structure (inconsistent number of columns), so read
    # it row-wise
    text = readLines(file.path(dir_projectDef, "ampliconsProject.txt"))

    # select lines with sample data
    sample_lines = grep("^Sample", text, value=TRUE)
    sample_con = textConnection(sample_lines)
    if(length(sample_lines) > 0){
    	sample_tab = read.csv(sample_con, sep="\t", col.names=c("Def", "Sample",
            "annotation", "name"), header=FALSE, colClasses=c("character",
            "character", "character", "character"), na.strings="", stringsAsFactors=FALSE)
	samples = sample_tab$Sample
	sampleID = sample_tab$name
    	numSamples = nrow(sample_tab)
    }else{
	samples = NA
    	warning(paste("sample information missing in", file.path(dir_projectDef, "ampliconsProject.txt")))
    }
    close(sample_con)

    if(!any(is.na(samples))){

    	# select lines with ReadData
    	RD_lines = grep("^RD\t", text, value=TRUE)
    	RD_con = textConnection(RD_lines)
    	if(length(RD_lines) > 0)
    	    RD_tab = read.csv(RD_con, sep="\t", col.names=c("Def", "RD", "active",
            	"annotation", "currentPath", "name", "originalPath", "readDataGroup", "sequenceBlueprint"),
            	header=FALSE, colClasses=c("character", "character", "character",
            	"character", "character", "character", "character", "character", "character"), na.strings="", stringsAsFactors=FALSE)
    	close(RD_con)

    	# select lines with ReadDataGroup
    	RDG_lines=grep("^RDG", text, value=TRUE)
    	RDG_con = textConnection(RDG_lines)
    	if(length(RDG_lines) > 0)
    	    RDG_tab = read.csv(RDG_con, sep="\t", col.names=c("Def", "RDG", "annotation",
            	"name"), header=FALSE, colClasses=c("character", "character",
            	"character", "character"), na.strings="", stringsAsFactors=FALSE)
    	close(RDG_con)

    	# select lines with lane ptp and read group information

    	RDS_lines = grep("^RD_Samp_Amp", text, value=TRUE)
    	RDS_con = textConnection(RDS_lines)
    	if(length(RDS_lines) > 0 & length(RD_lines) > 0 & length(RDG_lines) > 0){
    	    RDS_tab = read.csv(RDS_con, sep="\t", col.names=c("Def", "RD_Samp_Amp",
            	"amplicon", "readData", "sample"), header=FALSE, colClasses=c(
            	"character", "character", "character", "character", "character"), na.strings="", stringsAsFactors=FALSE)

    	    # filter those samples, which have several different ReadData-entries		
	    RDS_tab = subset(RDS_tab, duplicated(RDS_tab[,4:5]) == FALSE)

	    RData = merge(merge(RDS_tab, RD_tab, by.x="readData", by.y="RD", all.x=TRUE), RDG_tab, by.x="readDataGroup", by.y="RDG", all.x=TRUE)
	    MUX_MID1_tab = data.frame(name=rep(NA, numSamples), sample=samples)
	    MUX_MID2_tab = data.frame(name=rep(NA, numSamples), sample=samples)

    	}else{
    	    MID_lines = grep("^MIDSeq", text, value=TRUE)
    	    MID_con = textConnection(MID_lines)
    	    MUX_lines = grep("^Mux_Mid12_Samp", text, value=TRUE)
    	    MUX_con = textConnection(MUX_lines)
    	    RDMUX_lines = grep("^RD_Mux\t", text, value=TRUE)
    	    RDMUX_con = textConnection(RDMUX_lines)

    	    if(length(MID_lines) > 0 & length(MUX_lines) > 0 & length(RDMUX_lines) > 0 & length(RD_lines) > 0 & length(RDG_lines) > 0){
    	    	MID_tab = read.csv(MID_con, sep="\t", col.names=c("Def", "MIDSeq", "annotation", "midGroup", "name", "sequence"), 
		    header=FALSE, colClasses=c("character", "character", "character", "character", "character", "character"), na.strings="", stringsAsFactors=FALSE)
    	    	MUX_tab = read.csv(MUX_con, sep="\t", col.names=c("Def", "Mux_Mid12_Samp", "mid1", "mid2", "mux", "sample"), 
		    header=FALSE, colClasses=c("character", "character", "character", "character", "character", "character"), na.strings="", stringsAsFactors=FALSE)
    	    	RDMUX_tab = read.csv(RDMUX_con, sep="\t", col.names=c("Def", "RD_Mux", "mux", "readData"), 
		    header=FALSE, colClasses=c("character", "character", "character", "character"), na.strings="", stringsAsFactors=FALSE)

	    	MUX_MID1_tab = subset(merge(MUX_tab, MID_tab, by.x="mid1", by.y="MIDSeq", all.x=TRUE), !is.na(name))
	    	MUX_MID2_tab = subset(merge(MUX_tab, MID_tab, by.x="mid2", by.y="MIDSeq", all.x=TRUE), !is.na(name))
	    	RData = merge(merge(merge(MUX_tab, RDMUX_tab, by.x="mux", by.y="mux", all.x=TRUE, suffixes=c(".mux", ".rdmux")), RD_tab, by.x="readData", by.y="RD", 
		    all.x=TRUE), RDG_tab, by.x="readDataGroup", by.y="RDG", all.x=TRUE, suffixes=c(".rd", ".rdg"))

                RData = merge(merge(merge(MUX_tab, RDMUX_tab, by.x="mux", by.y="mux", all.x=TRUE, suffixes=c(".mux", ".rdmux")), RD_tab, by.x="readData", by.y="RD", 
		    all.x=TRUE), RDG_tab, by.x="readDataGroup", by.y="RDG", all.x=TRUE)


                
	    }else{
	    	warning(paste("Read data or MID entries missing in", file.path(dir_projectDef, "ampliconsProject.txt")))
	    	MUX_MID1_tab = data.frame(name=rep(NA, numSamples), sample=samples)
	    	MUX_MID2_tab = data.frame(name=rep(NA, numSamples), sample=samples)
	    	RData = data.frame(sample=samples, currentPath=rep(NA, numSamples), 
		     annotation.x=rep(NA, numSamples), name.y=rep(NA, numSamples))
	    }
	    close(MID_con)
    	    close(MUX_con)
    	    close(RDMUX_con)
    	}
    	close(RDS_con)

    	sampleMatrix = matrix(NA, numSamples, 7)
    	i = 1
    	for(s in samples){
	    mid1 = paste(subset(MUX_MID1_tab, sample==s)$name,collapse=",") 
	    mid2 = paste(subset(MUX_MID2_tab, sample==s)$name, collapse=",")
	    path = unique(c(subset(RData, sample==s)$currentPath, subset(RData, sample==s)$currentPath))
	    if(!any(is.na(path))){
	        path = sapply(strsplit(path, split="\\."), function(x)x[1])
	        ptp = paste(substr(path, 1, nchar(path)-2), collapse=",")
	        lane = paste(substr(path, nchar(path)-1, nchar(path)), collapse=",")
	    }else{
	    	ptp = NA
	    	lane = NA
	    }
	    rdAnnot = paste(unique(subset(RData, sample==s)$annotation.rd), collapse=",")
	    rgName = paste(unique(subset(RData, sample==s)$name.rdg), collapse=",")
	    sampleMatrix[i, ] = c(s, mid1, mid2, ptp, lane, rgName, rdAnnot)
	    i = i+1
    	}
    	sampleData = data.frame(sampleMatrix, row.names= sample_tab$name, stringsAsFactors=FALSE)
    	colnames(sampleData) = c("SampleID", "MID1", "MID2", "PTP_AccNum", "Lane", "ReadGroup", "Annotation")

    }else
	sampleData = data.frame(SampleID=NA, MID1=NA, MID2=NA, PTP_AccNum=NA, Lane=NA, ReadGroup=NA, Annotation=NA)

    # no label desription needed
    sampleMeta=data.frame(row.names=colnames(sampleData),
    	labelDescription=c("-","-","-","-","-","-","-"))

    return(list(sampleData,sampleMeta))
})


# readVariants runs through the subdirectory ".../Amplicons/Results/Variants/"
# and reads feature- and assay data for each variant 
# output: list containing
# variantForwCounts, the number of occurences of the variants in a reference
# sequence (read forwards)
# totalForwCount, the number of forward reads
# variantRevCounts, the number of occurences of the variants in a reference
# sequence (read backwards)
# totalRevCount, the number of reverse reads
# featureData, for each variant the name, canonical pattern and reference
# sequence
# featureMeta, some further information about the featureData

setMethod("readVariants",
    signature(dir_variants="character", samples="data.frame",
        refSeqs="DNAStringSet"),
    function(dir_variants, samples, refSeqs){

        # read the featureData for the variants (its name, canonical pattern and
        # reference sequence):
        variantDefs=read.table(file=file.path(dir_variants, "currentVariantDefs.txt"), sep="\t",
            header=TRUE, colClasses=c("character", "character", "character",
            "character"), stringsAsFactors=FALSE)	
        # read position and sequence information for each variant
        num_variants=length(variantDefs$CVariant)
        referenceSeq=matrix(0, num_variants, 1)
        v_splits=strsplit(variantDefs$canonicalPattern, split="[\\)\\(,-]!?")
        variantSeq=sapply(v_splits, function(x)x[3])
        deletions=sapply(v_splits, function(x) x[1] == "d")
        start=as.numeric(sapply(v_splits,function(x)x[2]))
        end=start
        end[deletions]=as.numeric(variantSeq[deletions])
        del_length=end[deletions] - start[deletions] + 1
        del_length=sapply(del_length, function(x) paste(rep("-",x), collapse=""))
        
        if(length(del_length) > 0)
           variantSeq[deletions] = del_length
           
        # get the bases from the reference sequences corresponding to the
        # variants
        for(i in 1:num_variants) {
            v_split=v_splits[[i]]
            refSeq_start=as.numeric(v_split[2])
            if(v_split[1] == "s") {
                refSeq_end=as.numeric(v_split[2])
            } else {
                if(v_split[1] == "d") {
                    refSeq_end=as.numeric(v_split[3])
	        } else {
                    stop("illegal canonical pattern")
                }
            }
            refSeq=refSeqs[names(refSeqs) == variantDefs$referenceSeq[i]]
            referenceSeq[i]=toString(subseq(refSeq[[1]], refSeq_start,
                refSeq_end))
        }
        # build name of the variant
        name=paste(start, ":", referenceSeq, "/", variantSeq, sep="")
        name[deletions]=paste(start[deletions], "-", end[deletions], ":DEL(",
            end[deletions] - start[deletions] + 1, ")", sep="")

        featureData=data.frame(row.names=variantDefs$CVariant, name=name,
            canonicalPattern=variantDefs$canonicalPattern,
            referenceSeqID=variantDefs$referenceSeq,
            start=start, 
            end=end,
            variantBase=variantSeq,
            referenceBases=referenceSeq,
            stringsAsFactors=FALSE)

        # construct meta data for the features
	# no label description needed
        featureMeta=data.frame(row.names=colnames(featureData),
            labelDescription=c("-","-","-","-","-","-","-"))

        # two separate matrices for forward and reverse reads each
        sampleNames=as.character(rownames(samples))
        num_samples=length(sampleNames)
        init = matrix(0, length(variantDefs$CVariant), num_samples)
	rownames(init) = variantDefs$CVariant
        colnames(init) = sampleNames
        variantForwCount = init
        totalForwCount = init
        variantRevCount = init
        totalRevCount = init

        # read all detections for all samples
        for(i in 1:num_samples){
            s_id=samples$SampleID[i]
            s_name=sampleNames[i]
            if(file.exists(file.path(dir_variants, s_id))){
                detections = dir(file.path(dir_variants, s_id), pattern=".txt$", ignore.case=FALSE)
                for(d in detections){
                    det = read.table(file.path(dir_variants, s_id, d), sep="\t",
                    	header=TRUE, 
			col.names=c("Def", "Det_Samp_Var", "canonicalVariant", 
			"forwardCount", "forwardDepth", "readType", "reverseCount", 
			"reverseDepth", "sample "),
			colClasses=c("character", "character",
                        "character", "numeric", "numeric", "numeric",
                        "numeric", "numeric", "character"),
			stringsAsFactors=FALSE)
		    if(nrow(det) > 0){
                    	det=subset(det, det$readType == 0)
                    	variantForwCount[det$canonicalVariant, s_name]=
                            det$forwardCount
                    	variantRevCount[det$canonicalVariant, s_name]=
                            det$reverseCount
                    	totalForwCount[det$canonicalVariant, s_name]=
                            det$forwardDepth
                    	totalRevCount[det$canonicalVariant, s_name]=
                            det$reverseDepth
		    }
                }
            } else {
                warning(paste("no variant data for ",s_name))
            }
        }
        return(list(variantForwCount, totalForwCount,
            variantRevCount, totalRevCount, featureData,
            featureMeta))
    }
)

          
# readReferenceSequences runs through the file
# ".../Amplicons/ProjectDef/ampliconsProject.txt" and reads the reference
# sequences
# output: refSeqAligned, an AlignedRead data frame with the reference sequences
# and additional information (name, id, gene)

setMethod("readReferenceSequences",

    signature(dir_projectDef="character"),
    function(dir_projectDef){
        # file has no data frame structure (inconsistent number of columns), so
        # read it row-wise
        pfLines=readLines(file.path(dir_projectDef, "ampliconsProject.txt"))

        # select lines with reference sequences (beginning with "RefSeq")
        refSeq_lines=grep("^RefSeq\t", pfLines, value=TRUE)

        # referenceSeqs above are in data frame format, so convert them to tables
        refSeq_con=textConnection(refSeq_lines)
        refSeq_tab=read.csv(refSeq_con, sep="\t", col.names=c("Def", "RefSeq",
            "annotation", "name", "sequence", "startPos"), header=FALSE,               
            colClasses=c("character", "character", "character", "character",
                "character", "numeric"), stringsAsFactors=FALSE)

        # build DNAStringSet for the reference sequence slot
        refSeqDNASet=DNAStringSet(refSeq_tab$sequence)
        names(refSeqDNASet)=refSeq_tab$RefSeq
        refSeqNames=refSeq_tab$name
        refSeqGenes=sapply(refSeqNames, function(x) strsplit(x, "_")[[1]][1] )
        refSeqData=data.frame(row.names=refSeqNames, name=refSeqNames,
            refSeqID=refSeq_tab$RefSeq, gene=refSeqGenes, stringsAsFactors=FALSE)
        refSeqAligned = AlignedRead(sread=refSeqDNASet, 
            id=BStringSet(refSeq_tab$RefSeq),
            position=as.integer(rep(1, length(refSeq_tab$RefSeq))),
            alignData=AlignedDataFrame(data=refSeqData,
                metadata=data.frame(row.names=colnames(refSeqData))))
        close(refSeq_con)
        return(refSeqAligned)
    }
)


# readAmplicons runs through the file
# ".../Amplicons/ProjectDef/ampliconsProject.txt" and reads feature- and assay
# data for each amplicon
# output: list containing
# numReadsForward, number of forward reads for a given amplicon
# numReadsReverse, number of reverse reads for a given amplicon
# featureData, for each amplicon the name, reference sequence, primer and
# end/start
# featureMeta, some useful information about the featureData

setMethod("readAmplicons",
    signature(dir_projectDef="character", dir_align="character", samples="data.frame"),
    function(dir_projectDef, dir_align, samples){

        # read the featureData for the amplicons (its name, reference sequence,
        # primer and end/start):
        # file has no data frame structure (inconsistent number of columns), so
        # read it row-wise
        pfLines=readLines(file.path(dir_projectDef, "ampliconsProject.txt"))

        # select lines with amplicons/primers (beginning with "Amp") 
        primer_lines=grep("^Amp\t", pfLines, value=TRUE)

        # amplicons above are in data frame format, so convert them to tables
        primer_con=textConnection(primer_lines)
        primer_tab=read.csv(		
            primer_con, sep="\t",
            col.names=c("Def", "Amp", "annotation", "name", "primer1", "primer2",
                "referenceSeq", "targetEnd", "targetStart"),
            header=FALSE,
            colClasses=c("character", "character", "character", "character",
                "character", "character", "character", "numeric", "numeric"),
	    stringsAsFactors=FALSE)
        # filter the desired feature data
        featureData=data.frame(row.names=primer_tab$name,
            ampID=primer_tab$Amp, 
            primer1=primer_tab$primer1,
            primer2=primer_tab$primer2,
            referenceSeqID=primer_tab$referenceSeq,
            targetEnd=primer_tab$targetEnd,
            targetStart=primer_tab$targetStart,
            stringsAsFactors=FALSE)

        # construct meta data for the features
        featureMeta=data.frame(row.names=colnames(featureData),
            labelDescription=c("-", "-", "-", "-", "-", "-"))

        close(primer_con)

        # count the amplicons for each sample
        refSeqsInd=grep("^RefSeq", pfLines)
        refSeqSplit=strsplit(pfLines[refSeqsInd], "\t")
        refSeqs=data.frame(
            RefSeqID=sapply(refSeqSplit, function(x) {x[2]}),
            RefSeqName=sapply(refSeqSplit, function(x) {x[4]}))

	## initialize assay data consisting of forward and reverse coverage
	sampleNames = rownames(samples)
	ampNames = rownames(featureData)
	init = matrix(0, length(sampleNames), length(ampNames))
	init = matrix(0, length(ampNames), length(sampleNames))
	rownames(init) = ampNames
	colnames(init) = sampleNames
        numReadsForward = init
        numReadsReverse = init

	## create named vector for fast connection between an ampid and its reference
	ampid2ref = rownames(featureData)
	names(ampid2ref) = featureData$ampID

	## read assay data
        for (s in 1:nrow(samples)) {
	    sampleName = sampleNames[s]
            thisSampleDir=file.path(dir_align, samples$SampleID[s])
            for (r in refSeqs$RefSeqID) {
		
                thisRefDir=file.path(thisSampleDir, r)
                alignFile=file.path(thisRefDir, paste(samples$SampleID[s],
                    "_vs_", r, ".ampIds.txt", sep=""))

                if (file.exists(alignFile)) {
                    seqCounts = read.table(file=alignFile, sep="\t", col.names=c("ampId", "forwardCount", "reverseCount"), stringsAsFactors=FALSE)
		    numReadsForward[ampid2ref[seqCounts$ampId], sampleName] = seqCounts[,"forwardCount"]
		    numReadsReverse[ampid2ref[seqCounts$ampId], sampleName] = seqCounts[,"reverseCount"]
                } 
            }
	}
        return(list(numReadsForward,numReadsReverse,featureData,featureMeta))	
    }
)


#################################################################################
########## ENSEMBL METHODS ######################################################
#################################################################################

setMethod("readAmpliconPositions",
    signature(object="AVASet", dnaSet="DNAStringSet", genes="character"), 
    function (object, dnaSet, genes = "", dataset = "hsapiens_gene_ensembl") {
        upDownStream = max(width(dnaSet))
        ensembl = useMart("ensembl", dataset = dataset)
        geneInfo = getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", 
            "start_position", "end_position", "strand", "chromosome_name"), 
            filter = "hgnc_symbol", values = unique(genes), mart = ensembl)

        geneSeqs = getSequence(id=geneInfo$ensembl_gene_id, type="ensembl_gene_id", 
            seqType="gene_exon_intron", mart=ensembl)
    
        upSeqs = getSequence(id=geneInfo$ensembl_gene_id, type="ensembl_gene_id", 
            seqType="gene_exon_intron", upstream=upDownStream, mart=ensembl)
        upSeqs = upSeqs[match(geneSeqs$ensembl_gene_id, upSeqs$ensembl_gene_id), ]
        upSeqs$gene_exon_intron = gsub(geneSeqs$gene_exon_intron, "", upSeqs$gene_exon_intron, fixed=TRUE)
    
        downSeqs = getSequence(id=geneInfo$ensembl_gene_id, type="ensembl_gene_id", 
            seqType="gene_exon_intron", downstream=upDownStream, mart=ensembl)    
        downSeqs = downSeqs[match(geneSeqs$ensembl_gene_id, downSeqs$ensembl_gene_id), ]
        downSeqs$gene_exon_intron = gsub(geneSeqs$gene_exon_intron, "", downSeqs$gene_exon_intron, fixed=TRUE)
        
        gSeqs = DNAStringSet(paste(upSeqs$gene_exon_intron,
            geneSeqs$gene_exon_intron, downSeqs$gene_exon_intron, sep=""))
        names(gSeqs) = geneSeqs$ensembl_gene_id

        numSeqs = length(dnaSet)

	numAllSeqs = length(sread(referenceSequences(object)))
        pos = rep(1, numAllSeqs)
        chr = rep(NA, numAllSeqs)
        strand = factor(rep(NA, numAllSeqs), levels=c("+", "-"))
        ensemblId = rep(NA, numAllSeqs)
	names(pos) = names(sread(referenceSequences(object)))
	names(chr) = names(sread(referenceSequences(object)))
	names(strand) = names(sread(referenceSequences(object)))
	names(ensemblId) = names(sread(referenceSequences(object)))
        for (i in 1:numSeqs) {
  	    refName = names(dnaSet)[i]
            amp = dnaSet[[i]]
            ind = geneInfo$hgnc_symbol == genes[i]
            gene = gSeqs[[geneInfo[ind, "ensembl_gene_id"]]]
        
            m = matchPattern(amp, gene)
            if (length(m) != 1) {
                stop(paste("Could not match amplicon:", names(dnaSet)[i]))
            }
            pos[refName] = geneInfo$start_position[ind] - upDownStream + start(m) - 1
            chr[refName] = geneInfo$chromosome_name[ind]
            ensemblId[refName] = geneInfo$ensembl_gene_id[ind]
            if (geneInfo$strand[ind] == 1) {
                strand[refName] = "+"
            } else {
                strand[refName] = "-"
            }
        }
        newAlignedDataFrame = data.frame(pData(alignData(referenceSequences(object))), 
            ensemblId = ensemblId)
        newAlignedData = AlignedRead(sread = sread(referenceSequences(object)), 
            id=id(referenceSequences(object)), chromosome=chr, strand=strand,
            position = as.integer(pos), alignData = AlignedDataFrame(data = newAlignedDataFrame, 
            metadata = data.frame(row.names = colnames(newAlignedDataFrame))))
        referenceSequences(object) = newAlignedData
        return(object)
    }
)


#################################################################################
########## DISPLAY, GETTER AND SETTER METHODS ###################################
#################################################################################


# this version of show adds some lines about the amplicons to the usual output
# from the eSet version

setMethod("show",
    signature(object="AVASet"),
    function(object){
        if(sum(variantFilterPerc(object)) == 0) {
            message("Variants: ")
        } else {
            message("Variants: ", 
		paste("(", names(variantFilterPerc(object))," filter = ", variantFilterPerc(object),")", sep=""))
        }
        callNextMethod(object)
        message("\nAmplicons: ")
        message("assayDataAmp:", paste(nrow(assayDataAmp(object)$forwCount), "features, ",
            ncol(assayDataAmp(object)$forwCount), "samples"))
        message("  element names:", assayDataElementNames(assayDataAmp(object)))  
        message("featureDataAmp: ")
        show(object@featureDataAmp)
        message("\nReference sequences: ")
        show(referenceSequences(object))
    }
)

# these versions of featureData and fData return the filtered featureData of the AVASet if the
# filter is set to a value > 0; otherwise it calls the original Biobase version

setMethod("fData",
    signature(object="AVASet"),
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
    signature(object="AVASet"),
    function(object){
        if(sum(variantFilterPerc(object)) > 0) {
            #filter the featureData
            return ( subset( object@featureData, rownames(object@featureData@data) %in% variantFilter(object) ) )
	} else {
            callNextMethod(object)
	}
    }
)


# this version of assayData returns the filtered assayData of the AVASet if the
# filter is set to a value > 0; otherwise it calls the original Biobase version

setMethod("assayData",
    signature(object="AVASet"),
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

# extend the Biobase version of the "["-method to the addditional amplicon slot

setMethod("[",
    signature("AVASet"),
    function(x, i, j, ..., drop=FALSE){
        # change the assayData and featureData of the additional amplicon slot
        if(missing(drop)) {
            drop=FALSE
        }
        if(missing(i) && missing(j)) {
            if(length(list(...)) != 0) { 
                stop("specify genes or samples to subset; use '", substitute(x),
                    "$", names(list(...))[[1]],
                    "' to access phenoData variables")
            }
            return(x)
        }

        # call original BioBase method to change variant data
        # change amplicon assayData manually
        orig=assayDataAmp(x)
        if(missing(i)) {
            x=callNextMethod(x, , j, drop=drop)
            orig=lapply(orig, function(obj) obj[ , j, ..., drop = drop])
        } else {
            if (missing(j)) {
                # call original BioBase method to change variant data
                x=callNextMethod(x, i, , drop=drop) 
            } else { 
                # call original BioBase method to change variant data
                x=callNextMethod(x,i,j,drop=drop)
                orig=lapply(orig, function(obj) obj[ , j, ..., drop = drop])
            }
        }
        assayDataAmp(x)=orig
        return(x)
    }
)

## provide another subset method for samples, variants AND amplicons ([] only works for
## samples and variants)
setMethod("subset",
    signature("AVASet"),
          
    function(x, subset, dimension="amplicons"){

      object = x
      if(!missing(subset)){

          if(dimension=="amplicons"){

            ## catch some bad input
            if(is.vector(subset, mode="character")){
             	if(!all(subset %in% rownames(fDataAmp(object)))){
                   stop("invalid argument: given amplicons not found")
             	}
            }else{
                if(is.vector(subset, mode="numeric")){
                    if(!all(subset %in% 1:nrow(fDataAmp(object)))){
                        stop("invalid argument: given indices out of range")
                    }
                }else{
                    if(is.vector(subset, mode="logical")){
                        if(length(subset) != nrow(fDataAmp(object))){
                            stop("invalid argument: given subset is of wrong length")
                        }
                    }else
                        stop("invalid argument: given subset of incorrect type")
                }
            }                  

            ## removing amplicons might affect some variants and samples:

            amps = fDataAmp(object)
            vars = fData(object)
            remainingAmps = amps[subset, ]

            ## remove variants from feature data
            remainingVars = vector(mode="character")
            i = 1
            ## Todo: substitute for loop by apply command
            for(v in rownames(vars)){
              var = vars[v, ]
              amps_var = subset(remainingAmps, referenceSeqID == var$referenceSeq)
              ## in case of one reference for more than one amplicon: test if variant lies within amplicon range
              if(any((var$start >= amps_var$targetStart) & (var$end <= amps_var$targetEnd))){
                remainingVars[i] = v
                i = i + 1
              }
            }
            newFeatureDataAmp = new("AnnotatedDataFrame", data=remainingAmps,
              varMetadata=varMetadata(featureDataAmp(object)))
            newFeatureData = new("AnnotatedDataFrame", data=fData(object)[remainingVars, ],
              varMetadata=varMetadata(featureData(object)))

            ## remove samples/amplicons/variants from assayData
            remainingAssayAmp1 = assayDataAmp(object)$forwCount[rownames(remainingAmps), , drop=FALSE]
            remainingAssayAmp2 = assayDataAmp(object)$revCount[rownames(remainingAmps), , drop=FALSE]
            idx1 = !apply(remainingAssayAmp1, 2, function(col) all(col == 0))
            idx2 = !apply(remainingAssayAmp2, 2, function(col) all(col == 0))
            idx = idx1 | idx2
            remainingSamples = colnames(remainingAssayAmp1)

            newAssayDataAmp=assayDataNew("list",
              forwCount=remainingAssayAmp1[, idx, drop=FALSE],
              revCount=remainingAssayAmp2[, idx, drop=FALSE])
            newAssayData=assayDataNew("list",
              variantForwCount=assayData(object)$variantForwCount[remainingVars, remainingSamples, drop=FALSE],
              totalForwCount=assayData(object)$totalForwCount[remainingVars, remainingSamples, drop=FALSE],
              variantRevCount=assayData(object)$variantRevCount[remainingVars, remainingSamples, drop=FALSE],
              totalRevCount=assayData(object)$totalRevCount[remainingVars, remainingSamples, drop=FALSE])

            ##remove samples from pheno data
            newPhenoData=new("AnnotatedDataFrame", data=pData(object)[remainingSamples, ],
              varMetadata=varMetadata(phenoData(object)))

            ## ubdate object
            assayData(object) = newAssayData
            assayDataAmp(object) = newAssayDataAmp
            featureData(object) = newFeatureData
            featureDataAmp(object) = newFeatureDataAmp
            phenoData(object) = newPhenoData

            return(object)

          }else{
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
                stop("invalid argument: dimension can be either \"samples\", \"variants\" or \"amplicons\"")
              }
            }
          }
        }else{
          warning("no subset given; AVASet object left unchanged")
          return(object)
        }
    }
)

# fDataAmp/featureDataAmp returns the amplicon data from the featureData slot
# (same as fData for variants from the original Biobase eSet)

setMethod("fDataAmp",
    signature(object="AVASet"),
    function(object){
        return(pData(object@featureDataAmp))
    }
)

setMethod("featureDataAmp",
    signature(object="AVASet"),
    function(object){
        return(object@featureDataAmp)
    }
)

# assayDataAmp returns the amplicon data from the assayData slot (same as
# assayData for variants from the original Biobase eSet)

setMethod("assayDataAmp",
    signature(object="AVASet"),
    function(object){
        return(object@assayDataAmp)
    }
)

setMethod("referenceSequences",
    signature(object="AVASet"),
    function(object){
        return(object@referenceSequences)
    }
)

setMethod("dirs",
    signature(object="AVASet"),
    function(object){
        return(object@dirs)
    }
)

setReplaceMethod("featureDataAmp",
    signature(object="AVASet", value="AnnotatedDataFrame"),
    function(object, value){
        object@featureDataAmp=value
        return(object)
    }
)

setReplaceMethod("assayDataAmp",
    signature(object="AVASet", value="AssayData"),
    function(object, value){
        object@assayDataAmp=value
        return(object)
    }
)

setReplaceMethod("referenceSequences",
    signature(object="AVASet", value="AlignedRead"),
    function(object, value){
        object@referenceSequences=value
        return(object)
    }
)
