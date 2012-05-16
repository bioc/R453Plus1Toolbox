# private function doing the annotation - object must be a proper data frame
.annotateVariantBlock <- function(object, ensembl, ensemblSnp) {
#    dataset = "hsapiens_gene_ensembl"
#    snpDataset = "hsapiens_snp"

    # create final annotation data structure
    variants = list()
    for (v in row.names(object)) {
        variants[[v]] = list(genes=data.frame(),
            transcripts=data.frame(),
            exons=data.frame(),
            snps=data.frame())
    }


    # variants as IRanges
    variantRanges = RangedData(IRanges(start=object$start, end=object$end),
        ID=row.names(object), strand=object$strand,
        seqSur=object$seqSur, seqRef=object$seqRef,
        seqMut=object$seqMut, space=object$chromosome)


    # fetch all genomic information from ENSEMBL
    martRegions = paste(object$chromosome, object$start, object$end, sep=":")
#    ensembl = useMart("ensembl", dataset=dataset)
    bmAttributes = c(
        "ensembl_gene_id",
        "chromosome_name",
        "start_position",
        "end_position",
        "strand",
        "ensembl_transcript_id",
        "transcript_start",
        "transcript_end",
        "ensembl_exon_id",
        "exon_chrom_start",
        "exon_chrom_end",
        "rank",
        "5_utr_start",
        "5_utr_end",
        "3_utr_start",
        "3_utr_end",
        "cds_start",
        "cds_end",
        "cds_length"
    )
    cat("Fetching gene structure information from Ensembl ... ")
    bmInfo = getBM(attributes=bmAttributes, filters="chromosomal_region",
        values=martRegions, mart=ensembl)
    if (nrow(bmInfo) == 0) {
        cat("done.\n")
        return(variants)
    }
    bmInfoNames = getBM(attributes=c("ensembl_transcript_id",
        "external_gene_id", "external_transcript_id"),
        filters="ensembl_transcript_id",
        values=unique(bmInfo$ensembl_transcript_id),
        mart=ensembl)
    mInd = match(bmInfo$ensembl_transcript_id,
        bmInfoNames$ensembl_transcript_id)
    bmInfo$external_gene_id = bmInfoNames$external_gene_id[mInd]
    bmInfo$external_transcript_id = bmInfoNames$external_transcript_id[mInd]
    cat("done.\n")
  
  
    # annotate gene information
    cat("Processing affected genes ... ")
    cols = c("ensembl_gene_id", "start_position",
        "end_position", "chromosome_name", "strand", "external_gene_id")
    geneDF = bmInfo[!duplicated(bmInfo[,"ensembl_gene_id",]), cols]
    geneRanges = RangedData(IRanges(start=geneDF$start_position,
        end=geneDF$end_position),
        ensembl_gene_id=geneDF$ensembl_gene_id,
        strand=geneDF$strand,
        external_gene_id=geneDF$external_gene_id,
        space=geneDF$chromosome_name)
    matchList = findOverlaps(query=ranges(geneRanges),
        subject=ranges(variantRanges))
    for (chr in names(matchList)) {
        numHits = nrow(as.matrix(matchList[[chr]]))
        for (hit in setdiff(0:numHits, 0)) {
            varInd = as.matrix(matchList[[chr]])[hit, "subjectHits"]
            geneInd = as.matrix(matchList[[chr]])[hit, "queryHits"]
            variants[[variantRanges[chr][["ID"]][varInd]]]$genes =
            rbind(variants[[variantRanges[chr][["ID"]][varInd]]]$genes,
                as.data.frame(geneRanges[chr])[geneInd,])
        }
    }
    cat("done.\n")
  

    # annotate transcript information
    cat("Processing affected transcripts ... ")
    cols = c("ensembl_gene_id", "ensembl_transcript_id",
        "transcript_start", "transcript_end", "chromosome_name",
        "strand", "external_transcript_id")
    transcriptDF = bmInfo[!duplicated(bmInfo[,"ensembl_transcript_id",]), cols]
    transcriptRanges = RangedData(IRanges(start=transcriptDF$transcript_start,
        end=transcriptDF$transcript_end),
        ensembl_gene_id=transcriptDF$ensembl_gene_id,
        ensembl_transcript_id=transcriptDF$ensembl_transcript_id,
        external_transcript_id=transcriptDF$external_transcript_id,
        strand=transcriptDF$strand,
        space=transcriptDF$chromosome_name)
    matchList = findOverlaps(query=ranges(transcriptRanges),
        subject=ranges(variantRanges))
    for (chr in names(matchList)) {
        numHits = nrow(as.matrix(matchList[[chr]]))
        for (hit in setdiff(0:numHits, 0)) {
            varInd = as.matrix(matchList[[chr]])[hit, "subjectHits"]
            transInd = as.matrix(matchList[[chr]])[hit, "queryHits"]
            variants[[variantRanges[chr][["ID"]][varInd]]]$transcripts =
            rbind(variants[[variantRanges[chr][["ID"]][varInd]]]$transcripts,
                as.data.frame(transcriptRanges[chr])[transInd,])
        }
    }
    cat("done.\n")
  

    # annotate exon information
    cat("Processing affected exons including coding regions and UTRs ... ")
    cols = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id",
        "exon_chrom_start", "exon_chrom_end", "rank",
        "5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end",
        "cds_start", "cds_end", "cds_length", "chromosome_name", "strand")
    exonDF = bmInfo[!duplicated(
        bmInfo[,c("ensembl_transcript_id","ensembl_exon_id")]), cols]
    exonRanges = RangedData(
        IRanges(start=exonDF$exon_chrom_start, end=exonDF$exon_chrom_end),
        ensembl_gene_id=exonDF$ensembl_gene_id,
        ensembl_transcript_id=exonDF$ ensembl_transcript_id,
        ensembl_exon_id=exonDF$ensembl_exon_id,
        rank=exonDF$rank,
        X5_utr_start=exonDF[,"5_utr_start"],
        X5_utr_end=exonDF[,"5_utr_end"],
        X3_utr_start=exonDF[,"3_utr_start"],
        X3_utr_end=exonDF[,"3_utr_end"],
        cds_start=exonDF$cds_start,
        cds_end=exonDF$cds_end,
        cds_length=exonDF$cds_length,
        strand=exonDF$strand,
        space=exonDF$chromosome_name)
    matchList = findOverlaps(query=ranges(exonRanges),
        subject=ranges(variantRanges))
    for (chr in names(matchList)) {
        numHits = nrow(as.matrix(matchList[[chr]]))
        for (hit in setdiff(0:numHits, 0)) {
            varInd = as.matrix(matchList[[chr]])[hit, "subjectHits"]
            exonInd = as.matrix(matchList[[chr]])[hit, "queryHits"]
            tmpExon = as.data.frame(exonRanges[chr])[exonInd,]
            tmpVar = as.data.frame(variantRanges[chr])[varInd,]

            cds = .codingRegion(tmpExon)
      
            # Is mutation in 5'-UTR?
            if (is.na(tmpExon$X5_utr_start)) {
                tmpExon$utr_5 = FALSE
            } else {
                mm = as.matrix(findOverlaps(
                    IRanges(tmpExon$X5_utr_start, tmpExon$X5_utr_end),
                    IRanges(tmpVar$start, tmpVar$end)))
                tmpExon$utr_5 = nrow(mm) == 1
            }
      
            # Is mutation in 3'-UTR?
            if (is.na(tmpExon$X3_utr_start)) {
                tmpExon$utr_3 = FALSE
            } else {
                mm = as.matrix(findOverlaps(
                    IRanges(tmpExon$X3_utr_start, tmpExon$X3_utr_end),
                    IRanges(tmpVar$start, tmpVar$end)))
                tmpExon$utr_3 = nrow(mm) == 1
            }
      
            # Is mutation in coding region?
            if (is.na(cds[1])) {
                tmpExon$coding = FALSE
            } else {
                mm = as.matrix(findOverlaps(
                    IRanges(cds[1], cds[2]),
                    IRanges(tmpVar$start, tmpVar$end)))
                tmpExon$coding = nrow(mm) == 1
            }
      
            # If mutation is in coding region, we compute affected codons
            if (tmpExon$coding) {
                if (tmpExon$strand == 1) {
                    relCdsPosStart = tmpVar$start - cds[1] + tmpExon$cds_start
                    relCdsPosEnd = tmpVar$end - cds[1] + tmpExon$cds_start
                    seqSur = DNAString(as.character(tmpVar$seqSur))
                    seqMut = DNAString(as.character(tmpVar$seqMut))
                } else {
                    relCdsPosStart = cds[2] - tmpVar$end + tmpExon$cds_start
                    relCdsPosEnd = cds[2] - tmpVar$start + tmpExon$cds_start
                    seqSur = reverseComplement(
                        DNAString(as.character(tmpVar$seqSur)))
                    seqMut = reverseComplement(
                        DNAString(as.character(tmpVar$seqMut)))
                }
                if (tmpVar$strand == "-") {
                    seqSur = reverseComplement(seqSur)
                    seqMut = reverseComplement(seqMut)
                }

                # substitution
                if (substr(tmpVar$seqRef, 1, 1) != "-" &
                        substr(tmpVar$seqMut, 1, 1) != "-") {
                    tmpExon$numCodonStart = ceiling(relCdsPosStart / 3)
                    tmpExon$numCodonEnd = ceiling(relCdsPosEnd / 3)
                    relCodonPosStart = relCdsPosStart -
                        ((tmpExon$numCodonStart - 1) * 3)
                    relCodonPosEnd = relCdsPosEnd -
                        ((tmpExon$numCodonEnd - 1) * 3)

                    # two bases at each end
                    #tmpExon$codonRef = as.character(substr(seqSur,
                    #    4-relCodonPosStart, length(seqSur)+1-relCodonPosEnd))
                    # three bases at each end
                    tmpExon$codonRef = as.character(substr(seqSur,
                        5-relCodonPosStart, length(seqSur)-relCodonPosEnd))
                    tmpExon$codonMut = tmpExon$codonRef
                    substr(tmpExon$codonMut, relCodonPosStart,
                        nchar(tmpExon$codonMut)-(3-relCodonPosEnd)) =
                        as.character(seqMut)

                    tmpExon$AminoRef = as.character(translate(
                        DNAString(tmpExon$codonRef)))
                    tmpExon$AminoMut = as.character(translate(
                        DNAString(tmpExon$codonMut)))

                # Deletion
                } else if (substr(tmpVar$seqMut, 1, 1) == "-") {
                    tmpExon$numCodonStart = ceiling(relCdsPosStart / 3)
                    tmpExon$numCodonEnd = ceiling(relCdsPosEnd / 3)
                    relCodonPosStart = relCdsPosStart -
                        ((tmpExon$numCodonStart - 1) * 3)
                    relCodonPosEnd = relCdsPosEnd -
                        ((tmpExon$numCodonEnd - 1) * 3)
                
                    tmpExon$codonRef = as.character(substr(seqSur,
                        5-relCodonPosStart, length(seqSur)-relCodonPosEnd))
                    tmpExon$codonMut = tmpExon$codonRef
                    substr(tmpExon$codonMut, relCodonPosStart,
                        nchar(tmpExon$codonMut)-(3-relCodonPosEnd)) =
                        as.character(seqMut)
                    
                    tmpExon$AminoRef = as.character(translate(
                        DNAString(tmpExon$codonRef)))
                    codonMutSeq = gsub("-", "", tmpExon$codonMut)
                    if (nchar(codonMutSeq) == 0) {
                        tmpExon$AminoMut = ""
                    } else if (nchar(codonMutSeq) == 3) {
                        tmpExon$AminoMut = as.character(translate(
                            DNAString(codonMutSeq)))
                    } else if (nchar(codonMutSeq) > 3) { # 4 bases remain
                        tmpExon$AminoMut = paste(as.character(translate(
                            DNAString(substr(codonMutSeq, 1, 3)))),
                            "+ frame-shift")
                    } else if (nchar(codonMutSeq) == 1) { # 1 base remains
                        nextCodon = paste(codonMutSeq,
                            as.character(substr(seqSur,
                            length(seqSur)-relCodonPosEnd + 1,
                            length(seqSur)-relCodonPosEnd + 2)), sep="")
                        nextAS = as.character(translate(DNAString(nextCodon)))
                        if (nextAS == "*") {
                            tmpExon$AminoMut = nextAS
                        } else {
                            tmpExon$AminoMut = paste(nextAS, "+ frame-shift")
                        }
                    } else if (nchar(codonMutSeq) == 2) { # 2 base remain
                        nextCodon = paste(codonMutSeq,
                            as.character(substr(seqSur,
                            length(seqSur) - relCodonPosEnd + 1,
                            length(seqSur) - relCodonPosEnd + 1)), sep="")
                        nextAS = as.character(translate(DNAString(nextCodon)))
                        if (nextAS == "*") {
                            tmpExon$AminoMut = nextAS
                        } else {
                            tmpExon$AminoMut = paste(nextAS, "+ frame-shift")
                        }
                    }

                # Insertion
                } else if (substr(tmpVar$seqRef, 1, 1) == "-") {
                    tmpExon$numCodonStart = ceiling(relCdsPosStart / 3)
                    tmpExon$numCodonEnd = ceiling(relCdsPosEnd / 3)
                    relCodonPosStart = relCdsPosStart -
                        ((tmpExon$numCodonStart - 1) * 3)
                    relCodonPosEnd = relCdsPosEnd -
                        ((tmpExon$numCodonEnd - 1) * 3)
                    if (relCodonPosStart == 3) { # => relCodonPosEnd == 1
                        tmpExon$numCodonStart = tmpExon$numCodonStart + 1
                        lIns = length(seqMut)
                        fillCodon = c(0, 2, 1)[(lIns %% 3)+1]
                        tmpExon$codonRef = as.character(substr(seqSur,
                            start=4, stop=3+lIns+fillCodon))
                        tmpExon$codonMut = tmpExon$codonRef
                        substr(tmpExon$codonMut, 1, lIns) = as.character(seqMut)
                        
                    } else if (relCodonPosStart == 2) {
                        lIns = length(seqMut)
                        fillCodon = c(0, 2, 1)[((lIns+2) %% 3)+1]
                        tmpExon$codonRef = as.character(substr(seqSur,
                            start=2, stop=3+lIns+fillCodon))
                        tmpExon$codonMut = tmpExon$codonRef
                        substr(tmpExon$codonMut, 3, lIns+2) =
                            as.character(seqMut)
                        
                    } else if (relCodonPosStart == 1) {
                        lIns = length(seqMut)
                        fillCodon = c(0, 2, 1)[((lIns+1) %% 3)+1]
                        tmpExon$codonRef = as.character(substr(seqSur,
                            start=3, stop=3+lIns+fillCodon))
                        tmpExon$codonMut = tmpExon$codonRef
                        substr(tmpExon$codonMut, 2, lIns+1) =
                            as.character(seqMut)
                    }                    

                    tmpExon$AminoMut = as.character(translate(
                        DNAString(tmpExon$codonMut)))
                    codonRefSeq = gsub("-", "", tmpExon$codonRef)
                    if (nchar(codonRefSeq) == 0) {
                        tmpExon$AminoRef = ""
                    } else if (nchar(codonRefSeq) == 3) {
                        tmpExon$AminoRef = as.character(translate(
                            DNAString(codonRefSeq)))
                    } else if (nchar(codonRefSeq) > 3) { # 4 bases remain
                        tmpExon$AminoRef = paste(as.character(translate(
                            DNAString(substr(codonRefSeq, 1, 3)))),
                            "+ frame-shift")
                    } else { # 1 or 2 bases remain
                        tmpExon$AminoRef = "frame-shift"
                    }
                }

            # Mutation is not in coding region
            } else {
                tmpExon$numCodonStart = NA
                tmpExon$numCodonEnd = NA
                tmpExon$codonRef = NA
                tmpExon$codonMut = NA
                tmpExon$AminoRef = NA
                tmpExon$AminoMut = NA
            }
      
            variants[[variantRanges[chr][["ID"]][varInd]]]$exons =
                rbind(variants[[variantRanges[chr][["ID"]][varInd]]]$exons,
                tmpExon[, c("space", "start", "end", "width", "ensembl_gene_id",
                "ensembl_transcript_id", "ensembl_exon_id", "rank",
                "strand", "utr_5", "utr_3", "coding", "numCodonStart",
                "numCodonEnd", "codonRef", "codonMut", "AminoRef", "AminoMut")])
        }
    }
    cat("done.\n")


    # download SNP information
#    ensemblSnp = useMart("snp", dataset=snpDataset)
    bmAttributes = c(
        "refsnp_id",
        "chr_name",
        "chrom_start",
        "chrom_strand",
        "allele")
    cat("Fetching SNP information from Ensembl ... ")
    bmInfo = getBM(attributes=bmAttributes, filters="chromosomal_region",
        values=martRegions, mart=ensemblSnp)
    cat("done.\n")
#    bmInfo$A = sapply(strsplit(bmInfo$allele, "/"), function(x) {x[1]})
#    bmInfo$B = sapply(strsplit(bmInfo$allele, "/"), function(x) {x[2]})


    # annotate SNP informations
    cat("Processing known SNPs ... ")
    snpDF = bmInfo
    if (nrow(snpDF) > 0) {
        snpDF$A = sapply(strsplit(snpDF$allele, "/"), function(x) {x[1]})
        snpDF$B = sapply(strsplit(snpDF$allele, "/"), function(x) {x[2]})        
        snpRanges = RangedData(IRanges(start=snpDF$chrom_start,
            end=snpDF$chrom_start + pmax(nchar(snpDF$A), nchar(snpDF$B)) - 1),
            refsnp_id=snpDF$refsnp_id,
            chrom_strand=snpDF$chrom_strand,
            allele=snpDF$allele,
            A=snpDF$A,
            B=snpDF$B,
            space=snpDF$chr_name)
        matchList = findOverlaps(query=ranges(snpRanges),
            subject=ranges(variantRanges))
        for (chr in names(matchList)) {
            numHits = nrow(as.matrix(matchList[[chr]]))
            for (hit in setdiff(0:numHits, 0)) {
                varInd = as.matrix(matchList[[chr]])[hit, "subjectHits"]
                snpInd = as.matrix(matchList[[chr]])[hit, "queryHits"]
                tmpSNP = as.data.frame(snpRanges[chr])[snpInd,]
                tmpVar = as.data.frame(variantRanges[chr])[varInd,]

                # We cannot process entries of class "HGMD_MUTATION", because
                # no alleles are given.
                if (tmpSNP$allele != "HGMD_MUTATION") {

                    # test whether the known SNP is the variant
                
                    sA = as.character(tmpSNP$A)
                    vA = as.character(tmpVar$seqRef)
                    if (substr(vA, start=1, stop=1) == "-") {
                        vA = "-"
                    }
                    sB = as.character(tmpSNP$B)
                    vB = as.character(tmpVar$seqMut)
                    if (substr(vB, start=1, stop=1) == "-") {
                        vB = "-"
                    }
                    if ((tmpSNP$chrom_strand == 1 & tmpVar$strand == "-") |
                        (tmpSNP$chrom_strand == -1 & tmpVar$strand == "+")) {
                        sA = as.character(complement(DNAString(sA)))
                        sB = as.character(complement(DNAString(sB)))
                    }                
                    if ((sA == vA & sB == vB & tmpSNP$start == tmpVar$start &
                            tmpSNP$end == tmpVar$end) |
                            (sA == vA & sB == vB & tmpSNP$start == tmpVar$end &
                            tmpSNP$end == tmpVar$start) |
                            (sB == vA & sA == vB & tmpSNP$start == tmpVar$start &
                            tmpSNP$end == tmpVar$end) |
                            (sB == vA & sA == vB & tmpSNP$start == tmpVar$end &
                            tmpSNP$end == tmpVar$start)) {
                        tmpSNP$identical = TRUE
                    } else {
                        tmpSNP$identical = FALSE
                    }
                } else { # HGMD_MUTATIONs
                    tmpSNP$identical = FALSE
                }

                variants[[variantRanges[chr][["ID"]][varInd]]]$snps =
                    rbind(variants[[variantRanges[chr][["ID"]][varInd]]]$snps,
                        tmpSNP)
            }
        }
    }
    cat("done.\n")
  
    annoVars = new("AnnotatedVariants")
    annotatedVariants(annoVars) = variants
    return(annoVars)
}



.annotateVariants_AVASet <- function(object){
    refSeqs=referenceSequences(object)

    mInd=match(fData(object)$referenceSeqID, as.character(id(refSeqs)))

    variantData=data.frame(
        start=NA,
        end=NA,
        chromosome=chromosome(refSeqs)[mInd],
        strand=strand(refSeqs)[mInd],
        seqRef=fData(object)$referenceBases,
        seqSur=substr(sread(refSeqs)[mInd], fData(object)$start - 3,
            fData(object)$end + 3),
        seqMut=fData(object)$variantBase,
        row.names=row.names(fData(object)),
        stringsAsFactors=FALSE
    )

    RefStart = position(refSeqs)[mInd]
    RefEnd = position(refSeqs)[mInd] + width(refSeqs)[mInd] - 1
    ind = variantData$strand == "+"

    # plus strand    
    variantData$start[ind] = RefStart[ind] + fData(object)$start[ind] - 1
    variantData$end[ind] = RefStart[ind] + fData(object)$end[ind] - 1
    # minus strand
    variantData$start[!ind] = RefEnd[!ind] - fData(object)$end[!ind] + 1
    variantData$end[!ind] = RefEnd[!ind] - fData(object)$start[!ind] + 1
    
    return(annotateVariants(variantData))
}


.annotateVariants_MapperSet <- function(object, bsGenome){

    chrs = as.character(fData(object)$chr)
    starts = fData(object)$start
    ends = fData(object)$end

    ## determine surroundings of +/-3 based in reference sequence
    if (missing(bsGenome)) {
        library("BSgenome.Hsapiens.UCSC.hg19")
        bsGenome = Hsapiens
    }    
    surr = vector(mode="character", length=length(chrs))
    for(i in 1:length(chrs))
	surr[i] = toString(subseq(bsGenome[[paste("chr", chrs[i], sep="")]], starts[i] - 3, ends[i] + 3))

    ## prepare data frame with variant data
    variantData = data.frame(
        start = starts,
        end = ends,
        chromosome = chrs,
        strand = fData(object)$strand,
        seqRef = fData(object)$referenceBases,
        seqSur = surr,
        seqMut = fData(object)$variantBase,
        row.names=row.names(fData(object)),
        stringsAsFactors=FALSE
    )

    return(annotateVariants(variantData))
}


.annotateVariants_data.frame <- function(object) {

    blockSize = 40
    dataset = "hsapiens_gene_ensembl"
    snpDataset = "hsapiens_snp"

    ensembl = useMart("ensembl", dataset=dataset)
    ensemblSnp = useMart("snp", dataset=snpDataset)
    
    av = new("AnnotatedVariants")
    i = 1
    while (i <= nrow(object)) {
        j = min(i + blockSize - 1, nrow(object))
        message(paste("Annotating variants ", i, " to ", j, ".", sep=""))
        subDf = object[i:j, ]
        subAv = .annotateVariantBlock(subDf, ensembl, ensemblSnp)
        annotatedVariants(av) = c(annotatedVariants(av), annotatedVariants(subAv))              
        i = i + blockSize
    }
    return(av)
}


setMethod("annotateVariants", signature=signature(object="AVASet", bsGenome="missing"),
    .annotateVariants_AVASet)

setMethod("annotateVariants", signature=signature(object="MapperSet", bsGenome="missing"),
    .annotateVariants_MapperSet)

setMethod("annotateVariants", signature=signature(object="MapperSet", bsGenome="BSgenome"),
    .annotateVariants_MapperSet)

setMethod("annotateVariants", signature=signature(object="data.frame", bsGenome="missing"),
    .annotateVariants_data.frame)






# private (not exported) function
#
# returns the absolute position of the coding region c(start, end)
# c(NA, NA) means that the exon has no coding region
.codingRegion <- function(exon) {

    cdsStart = NA
    cdsEnd = NA

    # some transcripts are not protein coding at all - fetch these cases first
    if (is.na(exon$cds_end) & is.na(exon$X5_utr_end) & is.na(exon$X3_utr_end)) {
        return(c(cdsStart, cdsEnd))
    }

    # whole exon is coding
    if (is.na(exon$X5_utr_end) & is.na(exon$X3_utr_end)) {
        cdsStart = exon$start
        cdsEnd = exon$end

    # exon has 5'-UTR
    } else if (!is.na(exon$X5_utr_end) & is.na(exon$X3_utr_end)) {
        if (exon$strand == 1) {
            if (exon$X5_utr_end == exon$end) {
                cdsStart = NA
                cdsEnd = NA
            } else {
                cdsStart = exon$X5_utr_end + 1
                cdsEnd = exon$end
            }
        } else {
            if (exon$X5_utr_start == exon$start) {
                cdsStart = NA
                cdsEnd = NA
            } else {
                cdsStart = exon$start
                cdsEnd = exon$X5_utr_start - 1
            }
        }

    # exon has 3'-UTR
    } else if (is.na(exon$X5_utr_end) & !is.na(exon$X3_utr_end)) {
        if (exon$strand == 1) {
            if (exon$X3_utr_start == exon$start) {
                cdsStart = NA
                cdsEnd = NA
            } else {
                cdsStart = exon$start
                cdsEnd = exon$X3_utr_start - 1
            }
        } else {
            if (exon$X3_utr_end == exon$end) {
                cdsStart = NA
                cdsEnd = NA
        } else {
                cdsStart = exon$X3_utr_end + 1
                cdsEnd = exon$end
        }
    }

    # exon has 5'-UTR and 3'-UTR => coding region must be in between
    } else if (!is.na(exon$X5_utr_end) & !is.na(exon$X3_utr_end)) {
        if (exon$strand == 1) {
            cdsStart = exon$X5_utr_end + 1
            cdsEnd =exon$X3_utr_start - 1
        } else {
            cdsStart = exon$X3_utr_end + 1
            cdsEnd = exon$X5_utr_start - 1
        }
    }
  
    return(c(cdsStart, cdsEnd))
}
