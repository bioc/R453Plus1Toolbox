########################################################################################################################
#######################################   Helper functions   ###########################################################
########################################################################################################################

.catchBadInput <- function(object, blocks, transcripts, sampleCols, minMut, dir, title){
  
    if(length(blocks) > 0)
        if(!is.vector(blocks, mode="character"))
	    stop("invalid argument: please specify a character vector of block names for each variant")
    if(length(transcripts) > 0)
        if(!is.vector(transcripts, mode="character"))
	    stop("invalid argument: please specify a character vector of containing Ensembl transcript-IDs")
    if(!is.vector(sampleCols, mode="character")){
	stop("invalid argument: lease specify a character vector of column names of the sample data (phenoData) of your
		AVASet/MapperSet object")}
    if(!is.numeric(minMut) | (minMut < 0) | (minMut > 100))
	stop("invalid argument: please specify a numeric vale for minMut between 0 and 100")
    if(!is.character(dir))
	stop("invalid argument: please specify your target directory as a string")
    if(!is.character(title))
        stop("invalid argument: title has to be of type character")

}

## sets one of the buttons ("Reference", "Sample" or "Quality") as missing and one as disabled, level for links from deeper folders
.createNavigation <- function(missing="", disabled="", level=c("", "", "")) {
	
	buttonBack = "<FORM><INPUT TYPE=\"button\" VALUE=\"Back\" onClick=\"history.go(-1);return true;\"></FORM>"

	if(missing == "Reference") {
		buttonReference = ""
	} else {
		buttonReference = paste("<FORM METHOD=\"LINK\" ACTION=\"",
			level[1], "index.html\"><INPUT TYPE=\"submit\" VALUE=\"Variant report by reference\"", 
			ifelse(disabled == "Reference", " disabled", ""), "></FORM>", sep="")
	}

	if(missing == "Sample") {
		buttonSample = ""
	} else {
		buttonSample = paste("<FORM METHOD=\"LINK\" ACTION=\"",
			level[2], "index_samples.html\"><INPUT TYPE=\"submit\" VALUE=\"Variant report by samples\"",
			ifelse(disabled == "Sample", " disabled", ""), "></FORM>", sep="")
	}

	if(missing == "Quality") {
		buttonQuality = ""
	} else {
		buttonQuality = paste("<FORM METHOD=\"LINK\" ACTION=\"",
			level[3], "qualityReport/qualityReport.html\"> <INPUT TYPE=\"submit\" VALUE=\"Quality report\"",
			ifelse(disabled == "Quality", " disabled", ""), "></FORM>", sep="")
	}
    
	return(
		paste("<table border=\"0\">
        <tr>
        <td>", buttonBack, "</td>
        <td>", buttonReference, "</td>
		<td>", buttonSample, "</td>
        <td>", buttonQuality, "</td>
        </tr>
        </table>", sep="")
	)

}

########################################################################################################################
#######################################   HTML report init for GSMSets   ###############################################
########################################################################################################################

.initHtmlReportGSMSet <- function(object, annot, blocks=c(), transcripts=c(), sampleCols, minMut=3, dir="HTMLReport", title="Summary"){

    ## catch some bad input and missing arguments
    if(missing(sampleCols))
        sampleCols = colnames(pData(object))
    .catchBadInput(object, blocks, transcripts, sampleCols, minMut, dir, title)
    if(missing(annot)){
        warning("no variant information specified: loading data from ensembl
            for all variants")
        annot=annotateVariants(object)
    }

    annot = annotatedVariants(annot)
    if(!is.null(blocks) & length(blocks) != length(annot)){
	warning("argument of illegal size: blocks were ommitted")
	blocks = NULL
    }

    vars_ID = intersect(names(annot), rownames(fData(object)))
    vars_withoutAnnot = setdiff(rownames(fData(object)), names(annot))
    
    ## check if annotated variants and the feature data of the AVASet are consistent
    if(length(vars_ID) == 0)
      stop("Variant annotations are not consistent with the AVASet: variant names in annotations differ from the names in the AVASet")
    else{
    	object = object[vars_ID, ]
	annot = annot[vars_ID]
    }
    if(length(vars_withoutAnnot) > 0)
	warning(length(vars_withoutAnnot), " variants lack of annotation and were omitted: ", paste(vars_withoutAnnot, collapse=", "))

    ## sample data
    variantPercs = round(getVariantPercentages(object, direction="both") * 100, digits=2)
    variantPercsForward = round(getVariantPercentages(object, direction="forward") * 100, digits=2)
    variantPercsReverse = round(getVariantPercentages(object, direction="reverse") * 100, digits=2)

    numReads = assayData(object)$totalForwCount + assayData(object)$totalRevCount
    numReadsMut = assayData(object)$variantForwCount + assayData(object)$variantRevCount
    numReadsForward = assayData(object)$totalForwCount
    numReadsMutForward = assayData(object)$variantForwCount
    numReadsReverse = assayData(object)$totalRevCount
    numReadsMutReverse = assayData(object)$variantRevCount

    ## report the reference gene (with exon information (from the
    ## referenceSequences-slot of the object))
    references_short = rep(NA, length(vars_ID))
    references = rep(NA, length(vars_ID))

    ## report maximum variant percentage
    max_perc=apply(variantPercs, 1, max)[vars_ID]

    ## report the variant name
    vars = rownames(fData(object)[vars_ID, ])
    names(vars) = vars_ID

    ## report base in variant and reference sequence and its start and end
    refSeq_base = fData(object)[vars_ID, ]$referenceBases
    ## only show deletions, that have max 5 symbols
    refSeq_base[nchar(refSeq_base) > 5]  = "-"
    var_base = fData(object)[vars_ID,]$variantBase
    ## only show deletions, that have max 5 symbols
    var_base[nchar(var_base) > 5]  = "-"
    seq_start=  fData(object)[vars_ID, ]$start
    seq_end = fData(object)[vars_ID, ]$end

    ## report absolute variant position in genome
    chromosome = fData(object)[vars_ID, "chr"]
    abs_start = fData(object)[vars_ID, "start"]
    abs_end = fData(object)[vars_ID, "end"]

    ## report type of variant (deletion, insertion point)
    mutation_type = vector(mode="character")
    for(i in vars_ID){
 	## insertion
	if((substr(fData(object)[i, "referenceBases"],1,1) == "-") | 
	    (nchar(fData(object)[i, "referenceBases"]) < nchar(fData(object)[i, "variantBase"]))
	) 
	    mutation_type[i] = "insertion"
	## deletion
	else 
	    if((substr(fData(object)[i, "variantBase"],1,1) == "-") | 
	    	(nchar(fData(object)[i, "referenceBases"]) > nchar(fData(object)[i, "variantBase"]))
	    ) 
		mutation_type[i] = "deletion"
	## point mutation
	    else mutation_type[i] = "point"
    }
    mutation_type = mutation_type[vars_ID]

    ## report change of amino acid
    as_change=as.character(sapply(annot, function(x) {
        if (nrow(x$exon) == 0) {
            return(FALSE)
        } else {
            aaChanges=x$exons$coding &
                ((x$exons$AminoRef != x$exons$AminoMut) |
                is.na(x$exons$AminoRef))
            return(sum(aaChanges, na.rm=TRUE) > 0)
        }
    }))
    as_change=sapply(as_change, function(x) if(x == "TRUE") x="yes"
        else x="no")

    ## report if variant is well known
    snps=sapply(annot, function(x)
        if(ncol(x$snps) > 0) {
	    paste(
            	paste("<a href=\"http://www.ensembl.org/Homo_sapiens/Variation/",
                     "Summary?source=dbSNP;v=", x$snps$refsnp_id,
                     "\" target=\"_blank\">", x$snps$refsnp_id, "</a>", sep=""),
		collapse=",")
        } else {
          "-"
        }
    )

    ## titles for the transcript pages
    genes = sapply(annot, function(x) {return(as.character(x$genes$external_gene_id[1]))})

    transTitles = vector(mode="character")    
    transTitles[vars_ID] = paste("Gene: ", genes[vars_ID], " - Variant: ", vars[vars_ID])

    ## build html report
    .htmlReport(vars_ID, variantPercs, variantPercsForward, variantPercsReverse, numReads, numReadsMut, 
	numReadsForward, numReadsMutForward, numReadsReverse, numReadsMutReverse,  max_perc, references, references_short, 
	vars, refSeq_base, var_base, chromosome, abs_start, abs_end, mutation_type, as_change, snps, transTitles,
	object, annot, blocks, transcripts, sampleCols, minMut, dir, title)
}

########################################################################################################################
#######################################   HTML report init for AVASets   ##############################################
########################################################################################################################

.initHtmlReportAVASet <- function(object, annot, blocks=c(), transcripts=c(), sampleCols, minMut=3, dir="HTMLReport", title="Summary"){
  
    ## catch some bad input and missing arguments
    if(missing(sampleCols))
        sampleCols = colnames(pData(object))

    .catchBadInput(object, blocks, transcripts, sampleCols, minMut, dir, title)
    if(missing(annot)){
        message("No variant information specified: loading data from ensembl for all variants")
        annot=annotateVariants(object)
    }

    annot = annotatedVariants(annot)
    if(!is.null(blocks) & length(blocks) != length(annot)){
	warning("argument of illegal size: blocks were ommitted")
	blocks = NULL
    }

    vars_ID = intersect(names(annot), rownames(fData(object)))
    vars_withoutAnnot = setdiff(rownames(fData(object)), names(annot))
    
    ## check if annotated variants and the feature data of the AVASet are consistent
    if(length(vars_ID) == 0)
      stop("Variant annotations are not consistent with the AVASet: variant names in annotations differ from the names in the AVASet")
    else{
    	object = object[vars_ID, ]
	annot = annot[vars_ID]
    }
    if(length(vars_withoutAnnot) > 0)
	warning(length(vars_withoutAnnot), " variants lack of annotation and were omitted: ", paste(vars_withoutAnnot, collapse=", "))
    
    ## sample data
    variantPercs = round(getVariantPercentages(object, direction="both") * 100, digits=2)
    variantPercsForward = round(getVariantPercentages(object, direction="forward") * 100, digits=2)
    variantPercsReverse = round(getVariantPercentages(object, direction="reverse") * 100, digits=2)

    numReads = assayData(object)$totalForwCount + assayData(object)$totalRevCount
    numReadsMut = assayData(object)$variantForwCount + assayData(object)$variantRevCount
    numReadsForward = assayData(object)$totalForwCount
    numReadsMutForward = assayData(object)$variantForwCount
    numReadsReverse = assayData(object)$totalRevCount
    numReadsMutReverse = assayData(object)$variantRevCount

    ## report the reference gene (with exon information (from the
    ## referenceSequences-slot of the object))
    references_short = sapply(annot, function(x) x$genes$external_gene_id)
    references=sapply(fData(object)[vars_ID, ]$referenceSeqID,
        function(x) subset(pData(alignData(referenceSequences(object))),
            pData(alignData(referenceSequences(object)))$refSeqID == x)$name)
    names(references) = vars_ID

    ## report maximum variant percentage
    max_perc=apply(variantPercs, 1, max)[vars_ID]

    ## report the variant name
    vars=fData(object)[vars_ID, ]$name
    names(vars)=vars_ID

    ## report base in variant and reference sequence and its start and end
    refSeq_base=fData(object)[vars_ID, ]$referenceBases
    ## only show deletions, that have max 5 symbols
    refSeq_base[nchar(refSeq_base) > 5]  = "-"
    var_base = fData(object)[vars_ID,]$variantBase
    ## only show deletions, that have max 5 symbols
    var_base[nchar(var_base) > 5]  = "-"
    seq_start=fData(object)[vars_ID, ]$start
    seq_end=fData(object)[vars_ID, ]$end

    ## report absolute variant position in genome
    refSeqInd=match(fData(object)[vars_ID, ]$referenceSeq,
        as.character(id(referenceSequences(object))))
    chromosome=as.character(
        chromosome(referenceSequences(object)[refSeqInd]))
    posStrand = strand(referenceSequences(object))[refSeqInd] == "+"
    abs_start = numeric(length(chromosome))
    abs_end = numeric(length(chromosome))
    abs_start[posStrand] = position(referenceSequences(object)[refSeqInd])[posStrand] + seq_start[posStrand] - 1
    abs_end[posStrand] = position(referenceSequences(object)[refSeqInd])[posStrand] + seq_end[posStrand] - 1
    endRefSeqs = position(referenceSequences(object)[refSeqInd]) + width(sread(referenceSequences(object)[refSeqInd])) - 1
    abs_start[!posStrand] = endRefSeqs[!posStrand] - seq_end[!posStrand] + 1
    abs_end[!posStrand] = endRefSeqs[!posStrand] - seq_start[!posStrand] + 1
    rm(refSeqInd)

    ## report type of variant (deletion, insertion point)
    mutation_type=sapply(fData(object)[vars_ID, ]$canonicalPattern,
        function(x) {substr(x, 1, 1)})
    mutation_type[mutation_type == "s"]="point"
    mutation_type[mutation_type == "d"]="deletion"
    mutation_type[mutation_type == "i"]="insertion"

    ## report change of amino acid
    as_change=as.character(sapply(annot, function(x) {
        if (nrow(x$exons) == 0) {
            return(FALSE)
        } else {
            aaChanges=x$exons$coding &
                ((x$exons$AminoRef != x$exons$AminoMut) |
                is.na(x$exons$AminoRef))
            return(sum(aaChanges, na.rm=TRUE) > 0)
        }
    }))
    as_change=sapply(as_change, function(x) if(x == "TRUE") x="yes"
        else x="no")

    ## report if variant is well known
    snps=sapply(annot, function(x)
        if(ncol(x$snps) > 0) {
	    paste(
            	paste("<a href=\"http://www.ensembl.org/Homo_sapiens/Variation/",
                     "Summary?source=dbSNP;v=", x$snps$refsnp_id,
                     "\" target=\"_blank\">", x$snps$refsnp_id, "</a>", sep=""),
		collapse=",")
        } else {
          "-"
        }

    )

    ## titles for the transcript pages
    transTitles = vector(mode="character")
    transTitles[vars_ID] = paste("Reference: ", references[vars_ID], " - Variant: ", vars[vars_ID])

    ## build html report
    .htmlReport(vars_ID, variantPercs, variantPercsForward, variantPercsReverse, numReads, numReadsMut, 
	numReadsForward, numReadsMutForward, numReadsReverse, numReadsMutReverse,  max_perc, references, references_short, 
	vars, refSeq_base, var_base, chromosome, abs_start, abs_end, mutation_type, as_change, snps, transTitles,
	object, annot, blocks, transcripts, sampleCols, minMut, dir, title)
}


########################################################################################################################
#######################################   HTML REPORT   ################################################################
########################################################################################################################
## The report is structured into two (MapperSet) or three (AVASet) parts containing variant and quality 
## information:
## (1) The main page sums up given variant information like the name, type, reference gene, position.
##     Using the argument "blocks", the main page can be individually structured by assigning a block name to 
##     each variant. 
##     The main page can be further structured by samples. For a given AVASet object, every sample links to 
##     another short quality report showing only the amplicon coverage for this sample.
## (2) Every variant on the main page links to a page with further details about the affected genes and transcripts
##     (e.g. Ensembl gene-IDs, transcript-IDs, codon sequences, changes of amino acids (if coding)).
## (3) Only in case of AVASet object: A quality report shows the coverage of every amplicon in forward and/or 
##     reverse direction. Further plots display the coverage by MID and PTP (if this information is given in the 
##     pheno data of the object).
########################################################################################################################

.htmlReport <- function(vars_ID, variantPercs, variantPercsForward, variantPercsReverse, numReads, numReadsMut,
        numReadsForward, numReadsMutForward, numReadsReverse, numReadsMutReverse,  max_perc, references, references_short,
        vars, refSeq_base, var_base, chromosome, abs_start, abs_end, mutation_type, as_change, snps, transTitles,
        object, annot, blocks, transcripts, sampleCols, minMut, dir, title){

    message("Creating variant report ... ", appendLF = FALSE)

    ## create directories
    if(!(file.exists(dir)))
        dir.create(dir, recursive=TRUE)
    if(!(file.exists(file.path(dir,"trans"))))
        dir.create(file.path(dir,"trans"), recursive=TRUE)
    if(!(file.exists(file.path(dir,"qualityReport"))))
        dir.create(file.path(dir,"qualityReport"), recursive=TRUE)

    ## load style sheet
    cssFile=file.path(find.package("R2HTML"), "samples/R2HTML.css")
    s=file.copy(cssFile, file.path(dir, "R2HTML.css"), overwrite=TRUE)
    ## load java script for sorting tables
    jsFile=file.path(.path.package("R453Plus1Toolbox"),
        "javascript/sorttable.js")
    js=file.copy(jsFile, file.path(dir, "sorttable.js"), overwrite=TRUE)


    ## create a html file for every transcript
	mutation_samples = vector(mode="character", length=length(vars_ID))
	names(mutation_samples) = vars_ID

    aminos = getAminoAbbr()
    trans_links = rep("", length(vars_ID))
    names(trans_links) = vars_ID

    for(v in vars_ID){
        var=annot[[v]]
        num_trans=length(var$transcripts$ensembl_transcript_id)
	    if(num_trans == 0)
		  trans_links[v]="not available"
	    else{
            transHTML=HTMLInitFile(file.path(dir,"trans"),
                filename=paste("trans",v,sep=""), extension="html",
                Title="Transcript Summary", CSSFile="../R2HTML.css")
            cat("<script src=\"../sorttable.js\"></script>", file=transHTML, append=TRUE)
            HTML.title(transTitles[v], HR=1, file=transHTML)
            ## add buttons
			if(class(object) == "AVASet") {
				cat(.createNavigation(level=c("../", "../", "../")), file=transHTML, append=TRUE)
			} else {
				cat(.createNavigation(missing="Quality", level=c("../", "../", "../")), file=transHTML, append=TRUE)
			}

           	## for the links in the summary page
            trans_links[v]=paste("<a href=\"trans/trans", v,
                ".html\">transcripts</a>", sep="")

            ## special case: all variants on introns (no exon data)
            if(ncol(var$exons) == 0){
                introns = rep(TRUE, num_trans)
                report_trans = data.frame(
                    row.names= seq(1, num_trans),
                    GeneID=paste(
                        c(rep("<a href=\"http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", num_trans)),
                        var$transcripts$ensembl_gene_id,
                        rep("\" target=\"_blank\">",num_trans),
                        var$transcripts$ensembl_gene_id,
                        rep("</a>", num_trans), sep=""),
		    	    TranscriptID_tmp = var$transcripts$ensembl_transcript_id,
                    TranscriptID=paste(c(
                        rep("<a href=\"http://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=", num_trans)),
                        var$transcripts$ensembl_transcript_id,
                        rep("\" target=\"_blank\">",num_trans),
                        var$transcripts$ensembl_transcript_id,
                        rep("</a>", num_trans), sep=""),
                    Coding="",
                    UTR5="",
                    UTR3="",
                    ExonRank="",
                    Codon="",
                    RefCodonSeq="",
                    MutCodonSeq="",
                    RefAmino="",
                    RefAminoLong="",
                    MutAmino="",
                    MutAminoLong="",
                    stringsAsFactors=FALSE
                )       
            } else { ## usual case: some/all variants on exons 
                var=merge.data.frame(var$transcripts, var$exons,
                    by="ensembl_transcript_id", all=TRUE)
                    num_trans=length(var$ensembl_transcript_id) # update number of transcript in case of more exons than transcripts
                    var$space.y=factor(var$space.y, 
                        levels=c(levels(var$space.y),"")
                )
                var$ensembl_gene_id.y = factor(var$ensembl_gene_id.y,
                  	levels=c(levels(var$ensembl_gene_id.y), ""))
                var$ensembl_exon_id = factor(var$ensembl_exon_id,
                   	levels=c(levels(var$ensembl_exon_id), ""))
                var[is.na(var)]=""
                introns=var$rank == ""
                var$numCodon=paste(var$numCodonStart, var$numCodonEnd, sep="-")
                var$numCodon[var$numCodonStart == var$numCodonEnd]=
                    as.character(var$numCodonStart[var$numCodonStart == var$numCodonEnd])
                ################# Detect Deletions
                ind=nchar(var$AminoRef) == 1 & nchar(var$AminoMut) == 1
                refAminoLong = rep("", length(ind))
                mutAminoLong = rep("", length(ind))
                refAminoLong[ind] = aminos[var$AminoRef[ind]]
                mutAminoLong[ind] = aminos[var$AminoMut[ind]]
                #################
                report_trans = data.frame(
                   	row.names=seq(1,num_trans), 
                   	GeneID=paste(c(rep("<a href=\"http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", num_trans)),
                        var$ensembl_gene_id.x,
                        rep("\" target=\"_blank\">",num_trans),
                        var$ensembl_gene_id.x, rep("</a>", num_trans), sep=""),
		    	    TranscriptID_tmp = var$ensembl_transcript_id,
                    TranscriptID = paste(c(rep("<a href=\"http://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=", num_trans)), 
			        var$ensembl_transcript_id, rep("\" target=\"_blank\">",num_trans), 
			        var$ensembl_transcript_id, rep("</a>", num_trans), sep=""),
                    Coding=var$coding,
                    UTR5=var$utr_5,
                    UTR3=var$utr_3,
                    ExonRank=paste(c(rep("<a href=\"http://www.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=", num_trans)), 
			        var$ensembl_gene_id.x, rep(";t=" ,num_trans), 
			        var$ensembl_transcript_id,rep("\" target=\"_blank\">",num_trans), var$rank, rep("</a>", num_trans), sep=""),
                    Codon=var$numCodon,
                    RefCodonSeq=var$codonRef,
                    MutCodonSeq=var$codonMut,
                    RefAmino=var$AminoRef,
                    RefAminoLong=refAminoLong,
                    MutAmino=var$AminoMut,
                    MutAminoLong=mutAminoLong,
                    stringsAsFactors=FALSE
                )		
            }

            report_trans[introns,"ExonRank"]="Intron"

            colnames(report_trans)[5]="5\'UTR"
            colnames(report_trans)[6]="3\'UTR"

    	    ## order transcript data frame according to a given list of transcript-IDs
    	    varTranscripts = intersect(transcripts, report_trans$TranscriptID_tmp)
            transOrder = which(report_trans$TranscriptID_tmp %in% varTranscripts)
            transOrder = c(transOrder, setdiff(1:num_trans, transOrder))
	    report_trans = report_trans[transOrder, ]
	    rownames(report_trans) = 1:num_trans
	    report_trans$TranscriptID_tmp = NULL

	    ## create HTML-file
            HTML(report_trans, file=transHTML, classtable="sortable")
			
            ## additional table for the sample data
            ## only use valid sample columns
            sampleCols = sampleCols[sampleCols %in% colnames(pData(object))]
            samples = data.frame(pData(object)[, sampleCols])

	    varFreq = t(variantPercs[v, ,drop =FALSE])
	    varFreqForward = t(variantPercsForward[v, ,drop =FALSE])
	    varFreqReverse = t(variantPercsReverse[v, ,drop =FALSE])

            samples$varFreq = varFreq
            samples$varFreqForward = varFreqForward
            samples$varFreqReverse = varFreqReverse	
            samples$Reads = t(numReads[v, , drop=FALSE])
            samples$ReadsMut = paste(t(numReadsMut[v, , drop=FALSE]), " (", varFreq, "%)", sep="") 
            samples$ReadsForward = t(numReadsForward[v, , drop=FALSE])
            samples$ReadsMutForward = paste(t(numReadsMutForward[v, , drop=FALSE]), " (", varFreqForward, "%)", sep="")
            samples$ReadsReverse = t(numReadsReverse[v, , drop=FALSE])
            samples$ReadsMutReverse = paste(t(numReadsMutReverse[v, , drop=FALSE]), " (", varFreqReverse, "%)", sep="")
            colnames(samples) = c(sampleCols, "VarFreq", "VarFreqForward", "VarFreqReverse",  
                "DepthTotal", "VarReads", "DepthForward", "VarReadsForward", "DepthReverse", "VarReadsReverse")
                
            ## for a AVASet object, every sample entry links to its quality report
            sampleNames = rownames(pData(object))
            sampleNums = 1:length(rownames(pData(object)))
            if(class(object) == "AVASet")
              rownames(samples) = paste("<a href=\"../qualityReport/", sampleNums, "/coverageSample_", sampleNums, ".html\">",sampleNames , "</a>", sep="")
            else
                rownames(samples) = sampleNames

            HTML(samples[, c(sampleCols, "DepthTotal", "VarReads", "DepthForward", "VarReadsForward", "DepthReverse", "VarReadsReverse")], 
                file=transHTML, classtable="sortable")
            HTMLEndFile(file=transHTML)

            ## compute column entries concerning samples having mutations in more than minMut% of the reads (for main table; see below)
            rownames(samples) = sampleNames
            mutation_samples[v] = toString( rownames(subset(samples, VarFreqForward >= minMut & VarFreqReverse >= minMut)) )
            if(mutation_samples[v] == "")
                mutation_samples[v] = "-"
	    }
    }


	## MAIN HTML 1 (no special structure or structured by given blocks)
        ## create main html file with gene overview and variant summary
	titleSummary = "Amplicon Summary"
        geneHTML = HTMLInitFile(dir, filename="index", extension="html",
            Title = titleSummary, CSSFile="R2HTML.css")
       	cat("<script src=\"sorttable.js\"></script>", file=geneHTML,
            append=TRUE)
        HTML.title(title, HR=1, file=geneHTML)
	## link to a page that structures the data by samples (see below)
	if(class(object) == "AVASet") {
		cat(.createNavigation(disabled="Reference"), file=geneHTML, append=TRUE)
	} else {
		cat(.createNavigation(missing="Quality", disabled="Reference"), file=geneHTML, append=TRUE)
	}
        ## create data frame to build a html table
        report_summary = data.frame( 
	        row.names = 1:length(references),
            Reference = references,
            Reference_tmp = as.character(references_short),
            Variant = vars,
            RefSeq = refSeq_base,
            VarSeq = var_base,
            AAChange = as_change,
            Mutation = mutation_type,
            KnownSNP = snps,
            Max = max_perc,
	        Min = mutation_samples,	    
            Position = paste("<a href=\"http://www.ensembl.org/Homo_sapiens/",
                "Location/View?db=core;r=", chromosome, ":",abs_start, "-",
                abs_end, "\" target=\"_blank\">", chromosome, ":", abs_start,
                "-", abs_end, "</a>", sep=""),
            Details = trans_links,
            stringsAsFactors = FALSE
        ) 	
	## remove column containing samples with a minumum amount of reads if the threshold (minMut) equals zero
	if(minMut == 0) report_summary$Min=NA

	colnames(report_summary) = c("Reference", "Reference_tmp", "Variant",
	    "RefSeq", "VarSeq", "AAChange", "Mutation", "KnownSNP", "Max%", 
	    paste("Min", minMut, "%", sep=""), "Position", "Details")

	## remove unneeded columns (columns consisting of NA-values, temporary columns)
	removedCols = which(is.na(report_summary), arr.ind=TRUE)
	removedCols = c(2, removedCols[, "col"])

        ## if there is no blocks-argument create just one big table
        ## otherwise create a single table for every given gene
        if(is.null(blocks)) {
            HTML(report_summary[-(removedCols)], file=geneHTML, classtable="sortable")
	## otherwise split the data frame into several blocks and create a html table for each block
        } else {
	    report_summary = report_summary[-(removedCols)]
            for(b in unique(blocks)){
		idx = which(blocks == b)
                HTML.title(b, file=geneHTML)
                HTML(report_summary[idx,],
                    file=geneHTML, classtable="sortable")
            }
        }
        HTMLEndFile(file=geneHTML)
    #save(report_summary, file="/tmp/rep_sum.RData")
    ## create reference overview table
    #for(ref in unique(

	## MAIN HTML 2 (structured by samples)
        ## create main html file with gene overview and variant summary
        geneHTML2 = HTMLInitFile(dir, filename="index_samples", extension="html",
            Title = titleSummary, CSSFile="R2HTML.css")
       	cat("<script src=\"sorttable.js\"></script>", file=geneHTML2,
            append=TRUE)
        HTML.title(title, HR=1, file=geneHTML2)
	## link to the main page and to a quality report (only for AVASet objects)
	if(class(object) == "AVASet") {
		cat(.createNavigation(disabled="Sample"), file=geneHTML2, append=TRUE)
	} else {
		cat(.createNavigation(missing="Quality", disabled="Sample"), file=geneHTML2, append=TRUE)
	}
        ## sample overview: add sample pData with anchors to each sample
        samplesWithLinks = pData(object)
        sampleNames = rownames(samplesWithLinks)
        sampleNums = 1:length(rownames(samplesWithLinks))
        sampleAnchors = paste("<a href=\"#anchor", sampleNums, " \">", sampleNames, "</a>", sep="")
        rownames(samplesWithLinks) = sampleAnchors
        ## add links to quality reports for a AVASet
        ## (extra column "amplicon coverage" for the pData table at the beginning and a link for each sample name)
	if(class(object) == "AVASet"){
            samplesWithLinks$QualityReport = paste("<a href=\"qualityReport/", sampleNums, "/coverageSample_", sampleNums, ".html\"> amplicon coverage </a>", sep="")
            reportSampleTitles = paste("<a href=\"qualityReport/", sampleNums, "/coverageSample_", sampleNums, ".html\">,<font color=\"#FFFFFF\">",sampleNames, "</font> </a>", sep="")
        }
        else
          reportSampleTitles = sampleNames
        names(reportSampleTitles) = sampleNames
        
        HTML(samplesWithLinks, file=geneHTML2, classtable="sortable")
        ## create data frame to build a html table for every sample
	for(sample in 1:length(sampleNames)){
	    ## only show those variants that are relevant for the sample (i.e. that meet the minMut requirement)
	    sampleVars = which(variantPercsForward[vars_ID, sample] > minMut & variantPercsReverse[vars_ID, sample] > minMut)
	    if(length(sampleVars) > 0){
                report_summary = data.frame( 
	    	    row.names = 1:length(references[sampleVars]),
            	    Reference = references[sampleVars],
            	    Reference_tmp = as.character(references_short)[sampleVars],
            	    Variant = vars[sampleVars],
            	    RefSeq = refSeq_base[sampleVars],
            	    VarSeq = var_base[sampleVars],
            	    AAChange = as_change[sampleVars],
            	    Mutation = mutation_type[sampleVars],
            	    KnownSNP = snps[sampleVars],
            	    Max = max_perc[sampleVars],
	    	    Min = mutation_samples[sampleVars],	    
            	    Position = paste("<a href=\"http://www.ensembl.org/Homo_sapiens/",
                        "Location/View?db=core;r=", chromosome, ":",abs_start, "-",
                        abs_end, "\" target=\"_blank\">", chromosome, ":", abs_start,
                        "-", abs_end, "</a>", sep="")[sampleVars],
            	    Details = trans_links[sampleVars],
            	    stringsAsFactors = FALSE
            	) 

	    	## remove column containing samples with a minumum amount of reads if the threshold (minMut) equals zero
	    	if(minMut == 0) report_summary$Min=NA

	    	colnames(report_summary) = c("Reference", "Reference_tmp", "Variant",
	            "RefSeq", "VarSeq", "AAChange", "Mutation", "KnownSNP", "Max%", 
	            paste("Min", minMut, "%", sep=""), "Position", "Details")

	    	## remove unneeded columns (columns consisting of NA-values, temporary columns)
	    	removedCols = which(is.na(report_summary), arr.ind=TRUE)
	    	removedCols = c(2, removedCols[, "col"])
	    }
            ## set anchor to jump to this sample from sample overview on top of the page
            cat(paste("<a name=\"anchor", sample, "\"></a>", sep=""), file=geneHTML2, append=TRUE)
	    ## set sample title (in case of AVASet: link to sample quality report)
            HTML.title(reportSampleTitles[sample], file=geneHTML2)

	    ## export data frame with variant information into the html file
            ## only if at least one variant passed the minMut filter
	    if(length(sampleVars) > 0)
            	HTML(report_summary[-(removedCols)], file=geneHTML2, classtable="sortable")

	}
        HTMLEndFile(file=geneHTML2)
        message("done")
        
	if(class(object) == "AVASet")
	    .qualityReport(object, dir)

    }

########################################################################################################################
#######################################   Quality report   #############################################################
########################################################################################################################

## quality report showing amplicon coverage by amplicons, mids and ptps
## only available for AVASets
.qualityReport <- function(object, dir){

    message("Creating quality report ... ", appendLF = FALSE)
    dirReport = "qualityReport"

    ## a script to hide elements like tables (e.g. when pressing a button)
    hideTableScript =
      "<script language=\"javascript\">
       function toggle_it(itemID){
       if ((document.getElementById(itemID).style.display == 'none'))
       {
         document.getElementById(itemID).style.display = 'inline';
       } else {
         document.getElementById(itemID).style.display = 'none';
       }
       }
       </script>"


    ## save amplicon assaydata as csv file
    write.csv(t(assayDataAmp(object)$forwCount), file=file.path(dir, dirReport, "ampCovForward.csv"))
    write.csv(t(assayDataAmp(object)$revCount), file=file.path(dir, dirReport, "ampCovReverse.csv"))
    write.csv(t(assayDataAmp(object)$forwCount) + t(assayDataAmp(object)$revCount), file=file.path(dir, dirReport, "ampCov.csv"))
              
    ## prepare pData for amplicon, filtered ptp and mid coverage plots and tables
    ## Amplicons
    pDataAMP = pData(object)
    sampleNames = rownames(pDataAMP)
	sampleNums = 1:length(rownames(pDataAMP)) 
    rownames(pDataAMP) = paste("<a href=\"", sampleNums, "/coverageSample_", sampleNums, ".html\">",sampleNames , "</a>", sep="")
    pDataAMP$QualityReport = paste("<a href=\"", sampleNums, "/coverageSample_", sampleNums, ".html\">amplicon coverage</a>", sep="")
	## MID
    if(!is.element("MID1", colnames(pData(object))) || all(is.na(pData(object)$MID1)) || all(pData(object)$MID1 == "NA")){
      mids = NA
      warning("No MIDs found. Left out amplicon coverage report for MIDs.")
    }
    else{
      ind = grep(",", pData(object)$MID1, fixed=TRUE)
      if (length(ind) != 0)
        pDataMID = pData(object[, -ind])
      else
        pDataMID = pData(object)
	  pDataMID$sampleNumber = 1:length(rownames(pDataMID))  ## retain sample numbers beyond ordering
      pDataMID = pDataMID[order(pDataMID$MID1), ]    ## sort samples by mid
      mids = unique(pDataMID$MID1)
      sampleNames = rownames(pDataMID)
	  sampleNums = pDataMID$sampleNumber
	  pDataMID$sampleNumber <- NULL
      rownames(pDataMID) = paste("<a href=\"", sampleNums, "/coverageSample_", sampleNums, ".html\">",sampleNames , "</a>", sep="")
      pDataMID$QualityReport = paste("<a href=\"", sampleNums, "/coverageSample_", sampleNums, ".html\">amplicon coverage</a>", sep="")
    }
    ## PTP
    if(!is.element("PTP_AccNum", colnames(pData(object))) || all(is.na(pData(object)$PTP_AccNum))  || all(pData(object)$PTP_AccNum == "NA") ){
      ptps = NA
      warning("No PTP information found. Left out amplicon coverage report for PTPs.")
    }
    else{
      ind = union(grep(",", pData(object)$PTP_AccNum, fixed=TRUE), grep(",", pData(object)$Lane, fixed=TRUE))
      if (length(ind) != 0)
        pDataPTP = pData(object[, -ind])
      else
        pDataPTP = pData(object)
      pDataPTP$sampleNumber = 1:length(rownames(pDataPTP))  ## retain sample numbers beyond ordering
      pDataPTP = pDataPTP[order(pDataPTP$PTP_AccNum), ]
      ptps = unique(paste(as.character(pDataPTP$PTP_AccNum), as.character(pDataPTP$Lane), sep=""))
      sampleNames = rownames(pDataPTP)
	  sampleNums = pDataPTP$sampleNumber
      pDataPTP$sampleNumber <- NULL
      rownames(pDataPTP) = paste("<a href=\"", sampleNums, "/coverageSample_", sampleNums, ".html\">",sampleNames , "</a>", sep="")
      pDataPTP$QualityReport = paste("<a href=\"", sampleNums, "/coverageSample_", sampleNums, ".html\">amplicon coverage</a>", sep="")
    }
    
    ## set size of output images (heuristic, depending on the number of amplicons, mids ans ptps)
    numAmps = nrow(fDataAmp(object))
    pngWidth = max(300, 30 * (numAmps + 5))
    pngHeight = 600
    pdfWidth = pngWidth / 72
    pdfHeight = pngHeight / 72
    if(!all(is.na(mids))){
      numMIDs = length(mids)
      pngWidthMID = max(300, 30 * (numMIDs + 5))
      pngHeightMID = 600
      pdfWidthMID = pngWidthMID / 72
      pdfHeightMID = pngHeightMID / 72
    }
    if(!all(is.na(ptps))){
      numPTPs = length(ptps)
      pngWidthPTP = max(300, 30 * (numPTPs + 5))
      pngHeightPTP = 600
      pdfWidthPTP = pngWidthPTP / 72
      pdfHeightPTP = pngHeightPTP / 72
    }

    ## Coverage in both directions (together)
    filenameCov1 = "coverageAmp1"
    ## PNG file
    png(filename=file.path(dir, dirReport, paste(filenameCov1, ".png", sep="")), width=pngWidth, height=pngHeight)
    plotAmpliconCoverage(avaSet=object, bothDirections=FALSE, type="amplicon")
    dev.off()
    ## PDF file
    pdf(file=file.path(dir, dirReport, paste(filenameCov1, ".pdf", sep="")), width=pdfWidth, height=pdfHeight)
    plotAmpliconCoverage(avaSet=object, bothDirections=FALSE, type="amplicon")
    dev.off()
    ## Coverage in both directions (separated)
    filenameCov2 = "coverageAmp2"
    ## PNG file
    png(filename=file.path(dir, dirReport, paste(filenameCov2, ".png", sep="")), width=pngWidth*2, height=pngHeight)
    plotAmpliconCoverage(avaSet=object, bothDirections=TRUE, type="amplicon")
    dev.off()
    ## PDF file
    pdf(file=file.path(dir, dirReport, paste(filenameCov2, ".pdf", sep="")), width=pdfWidth*2, height=pdfHeight)
    plotAmpliconCoverage(avaSet=object, bothDirections=TRUE, type="amplicon")
    dev.off()

    if(!all(is.na(mids))){
        ## Coverage in both directions by MID (together)
        filenameCov3 = "coverageMID1"
        ## PNG file
        png(filename=file.path(dir, dirReport, paste(filenameCov3, ".png", sep="")), width=pngWidthMID, height=pngHeightMID)
        plotAmpliconCoverage(avaSet=object, bothDirections=FALSE, type="mid")
        dev.off()
        ## PDF file
        pdf(file=file.path(dir, dirReport, paste(filenameCov3, ".pdf", sep="")), width=pdfWidthMID, height=pdfHeightMID)
        plotAmpliconCoverage(avaSet=object, bothDirections=FALSE, type="mid")
        dev.off()
        ## Coverage in both directions by MID (separated)
        filenameCov4 = "coverageMID2"
        ## PNG file
        png(filename=file.path(dir, dirReport, paste(filenameCov4, ".png", sep="")), width=pngWidthMID*2, height=pngHeightMID)
        plotAmpliconCoverage(avaSet=object, bothDirections=TRUE, type="mid")
        dev.off()
        ## PDF file
        pdf(file=file.path(dir, dirReport, paste(filenameCov4, ".pdf", sep="")), width=pdfWidthMID*2, height=pdfHeightMID)
        plotAmpliconCoverage(avaSet=object, bothDirections=TRUE, type="mid")
        dev.off()
    }

    if(!all(is.na(ptps))){
        ## Coverage in both directions by PTP (together)
        filenameCov5 = "coveragePTP1"
        ## PNG file
        png(filename=file.path(dir, dirReport, paste(filenameCov5, ".png", sep="")), width=pngWidthPTP, height=pngHeightPTP)
        plotAmpliconCoverage(avaSet=object, bothDirections=FALSE, type="ptp")
        dev.off()
        ## PDF file
        pdf(file=file.path(dir, dirReport, paste(filenameCov5, ".pdf", sep="")), width=pdfWidthPTP, height=pdfHeightPTP)
        plotAmpliconCoverage(avaSet=object, bothDirections=FALSE, type="ptp")
        dev.off()
        ## Coverage in both directions by PTP (separated)
        filenameCov6 = "coveragePTP2"
        ## PNG file
        png(filename=file.path(dir, dirReport, paste(filenameCov6, ".png", sep="")), width=pngWidthPTP*2, height=pngHeightPTP)
        plotAmpliconCoverage(avaSet=object, bothDirections=TRUE, type="ptp")
        dev.off()
        ## PDF file
        pdf(file=file.path(dir, dirReport, paste(filenameCov6, ".pdf", sep="")), width=pdfWidthPTP*2, height=pdfHeightPTP)
        plotAmpliconCoverage(avaSet=object, bothDirections=TRUE, type="ptp")
        dev.off()
    }


    ## coverage for each sample
    samples = rownames(pData(object))
    for(s in 1:length(samples)){
	dirSample = file.path(dirReport, s)
        if(!(file.exists(file.path(dir, dirSample))))
            dir.create(file.path(dir, dirSample), recursive=TRUE)
	## PNG file
    	png(filename=file.path(dir, dirSample, paste("coverageSample_", s, ".png", sep="")), width=pngWidth, height=pngHeight)
        plotAmpliconCoverage(avaSet=object[, s], bothDirections=TRUE)
	dev.off()
	## PDF file
    	pdf(file=file.path(dir, dirSample, paste("coverageSample_", s, ".pdf", sep="")), width=pdfWidth, height=pdfHeight)
        plotAmpliconCoverage(avaSet=object[, s], bothDirections=TRUE)
	dev.off()
    }

    ## HTML file for general amplicon coverage (all samples altogether)
    qualityHTML = HTMLInitFile(file.path(dir, dirReport), filename="qualityReport", extension="html",
        Title = "Quality report", CSSFile="../R2HTML.css")
    cat("<script src=\"../sorttable.js\"></script>", file=qualityHTML, append=TRUE)
    cat(hideTableScript, file=qualityHTML, append=TRUE)
    
    ## Amplicon Coverage
    HTML.title("Quality report", HR=1, file=qualityHTML)

    ## add navigation buttons
	cat(.createNavigation(disabled="Quality", level=c("../", "../", "../")), file=qualityHTML, append=TRUE)

	## add page internal navigation buttons
	if(!all(is.na(mids)))
      buttonCovMID =
        "<FORM METHOD=\"LINK\" ACTION=\"#anchorCovMID\">
	 <INPUT TYPE=\"submit\" VALUE=\"Amplicon coverage by MID\">
	 </FORM>"
    else
      buttonCovMID = ""
    if(!all(is.na(ptps)))
      buttonCovPTP =
        "<FORM METHOD=\"LINK\" ACTION=\"#anchorCovPTP\">
	 <INPUT TYPE=\"submit\" VALUE=\"Amplicon coverage by PTP\">
	 </FORM>"
    else
      buttonCovPTP = ""
	if(!all(is.na(mids)) | !all(is.na(ptps))) {
		cat(paste("Go to:<br />
		<table border=\"0\">
        <tr>
        <td>", buttonCovMID, "</td>
        <td>", buttonCovPTP, "</td>
        </tr>
        </table>", sep = ""), file=qualityHTML, append=TRUE)
	}

    ## add some links for download of assaydata
    cat("Download .csv files of:", file=qualityHTML, append=TRUE)
    cat("<br />", file=qualityHTML, append=TRUE)    
    cat("<a href=\"ampCovForward.csv\"> amplicon coverage (forward direction) </a>", file=qualityHTML, append=TRUE)
    cat("<br />", file=qualityHTML, append=TRUE)    
    cat("<a href=\"ampCovReverse.csv\"> amplicon coverage (reverse direction) </a>", file=qualityHTML, append=TRUE)
    cat("<br />", file=qualityHTML, append=TRUE)    
    cat("<a href=\"ampCov.csv\"> amplicon coverage (both directions) </a>", file=qualityHTML, append=TRUE)

    border = 0
    ## Amplicon coverage by amplicons
    HTML.title("Amplicon coverage",file=qualityHTML)
    ## insert image with link to a pdf version
    cat(paste("<a href=\"", filenameCov1, ".pdf\" target=\"_blank\">
	<img src=\"", filenameCov1, ".png\" width=", pngWidth, "height=", pngHeight, "border=", border, "alt=Coverage />
	</a>", sep=""), file=qualityHTML, append=TRUE)
    cat("<br />", file=qualityHTML, append=TRUE)
    ## insert image with link to a pdf version
    cat(paste("<a href=\"", filenameCov2, ".pdf\" target=\"_blank\">
	<img src=\"", filenameCov2, ".png\" width=", pngWidth*2, "height=", pngHeight, "border=", border, "alt=Coverage />
	</a>", sep=""), file=qualityHTML, append=TRUE)
    ## insert table with sample information
    ## (add possibility to hide the table; No better solution  yet for smuggling of table id and style in Border-value)
    buttonHideTab =
      "<FORM> <INPUT TYPE=\"button\" VALUE=\"Show/hide sample information\" onClick=\"toggle_it('sampleTabAmp')\"> </FORM>"
    cat(buttonHideTab, file=qualityHTML, append=TRUE)
    HTML(pDataAMP, file=qualityHTML, classtable="sortable", Border="1 id=\"sampleTabAmp\" style=\"display:none;\"")
    
    ## Amplicon Coverage by MIDs
    if(!all(is.na(mids))){
        cat("<a name=\"anchorCovMID\"></a>", file=qualityHTML, append=TRUE)
        HTML.title("Amplicon coverage by MID",file=qualityHTML)
        ## insert image with link to a pdf version
        cat(paste("<a href=\"", filenameCov3, ".pdf\" target=\"_blank\">
            <img src=\"", filenameCov3, ".png\" width=", pngWidthMID, "height=", pngHeightMID, "border=", border, "alt=Coverage />
            </a>", sep=""), file=qualityHTML, append=TRUE)
        cat("<br />", file=qualityHTML, append=TRUE)
        ## insert image with link to a pdf version
        cat(paste("<a href=\"", filenameCov4, ".pdf\" target=\"_blank\">
            <img src=\"", filenameCov4, ".png\" width=", pngWidthMID*2, "height=", pngHeightMID, "border=", border, "alt=Coverage />
            </a>", sep=""), file=qualityHTML, append=TRUE)
        cat("<br />", file=qualityHTML, append=TRUE)
        ## insert table with sample information
        ## (add possibility to hide the table; No better solution  yet for smuggling of table id and style in Border-value)
        buttonHideTab =
          "<FORM> <INPUT TYPE=\"button\" VALUE=\"Show/hide sample information\" onClick=\"toggle_it('sampleTabMID')\"> </FORM>"
        cat(buttonHideTab, file=qualityHTML, append=TRUE)
        HTML(pDataMID, file=qualityHTML, classtable="sortable", Border="1 id=\"sampleTabMID\" style=\"display:none;\"")
    }
    ## Amplicon Coverage by PTPs    
    if(!all(is.na(ptps))){
        cat("<a name=\"anchorCovPTP\"></a>", file=qualityHTML, append=TRUE)
        HTML.title("Amplicon coverage by PTP",file=qualityHTML)
        ## insert image with link to a pdf version
        cat(paste("<a href=\"", filenameCov5, ".pdf\" target=\"_blank\">
            <img src=\"", filenameCov5, ".png\" width=", pngWidthPTP, "height=", pngHeightPTP, "border=", border, "alt=Coverage />
            </a>", sep=""), file=qualityHTML, append=TRUE)
        cat("<br />", file=qualityHTML, append=TRUE)
            ## insert image with link to a pdf version
        cat(paste("<a href=\"", filenameCov6, ".pdf\" target=\"_blank\">
            <img src=\"", filenameCov6, ".png\" width=", pngWidthPTP*2, "height=", pngHeightPTP, "border=", border, "alt=Coverage />
            </a>", sep=""), file=qualityHTML, append=TRUE)
        cat("<br />", file=qualityHTML, append=TRUE)
        ## insert table with sample information
        ## (add possibility to hide the table; No better solution  yet for smuggling of table id and style in Border-value)
        buttonHideTab =
          "<FORM> <INPUT TYPE=\"button\" VALUE=\"Show/hide sample information\" onClick=\"toggle_it('sampleTabPTP')\"> </FORM>"
        cat(buttonHideTab, file=qualityHTML, append=TRUE)
        HTML(pDataPTP, file=qualityHTML, classtable="sortable", Border="1 id=\"sampleTabPTP\" style=\"display:none;\"")
    }
    HTMLEndFile(file=qualityHTML)

    ## Amplicon coverage for every sample (create html file for each sample)
    for(sample in 1:length(samples)){

		dirSample = file.path(dirReport, sample)
        sampleQualityHTML = HTMLInitFile(file.path(dir, dirSample), filename=paste("coverageSample_", sample, sep=""), 
          extension="html", Title = "Quality report", CSSFile="../../R2HTML.css")
        cat("<script src=\"../../sorttable.js\"></script>", file=sampleQualityHTML, append=TRUE)
    	HTML.title(paste("Amplicon coverage for sample", sample), HR=1, file=sampleQualityHTML)
    	## link to the main pages
		cat(.createNavigation(level=c("../../", "../../", "../../")), file=sampleQualityHTML, append=TRUE)

        ## insert image with link to a pdf version
        cat(paste("<a href=\"coverageSample_", sample, ".pdf\" target=\"_blank\">
            <img src=\"coverageSample_", sample, ".png\" width=", pngWidth, "height=", pngHeight, "border=", border, "alt=Coverage />
            </a>", sep=""), file=sampleQualityHTML, append=TRUE)
        HTMLEndFile(file=sampleQualityHTML)

    }
    message("done")

}

setMethod("htmlReport",
    signature=signature(object="AVASet"),
    .initHtmlReportAVASet)

setMethod("htmlReport",
    signature=c(object="MapperSet"),
    .initHtmlReportGSMSet)


