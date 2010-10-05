.alignShortReads_DNAStringSet <- function(object, bsGenome, seqNames, ensemblNotation) {

    if (is.null(names(object))) {
        names(object) = as.character(1:length(object))
    }

    pDictPlus = PDict(object, max.mismatch=NA,
        tb.start=1, tb.end=min(width(object)))
    pDictMinus = PDict(reverseComplement(object), max.mismatch=NA,
        tb.start=1, tb.end=min(width(object)))

    matches = data.frame(
        chr=as.character(rep(NA, length(object))),
        start=as.integer(rep(NA, length(object))),
        strand=as.character(rep(NA, length(object))),
        unique=rep(TRUE, length(object)),
        stringsAsFactors=FALSE)
    row.names(matches) = names(object)
  
    for (seqname in seqNames) {
        subject = bsGenome[[seqname]]
        index = list(
            "+"=matchPDict(pDictPlus, subject),
            "-"=matchPDict(pDictMinus, subject))

        for (s in c("+", "-")) {
            sindex = index[[s]]
            for (n in names(sindex)) {
                if (length(sindex[[n]]) == 1) {
                    if (matches[n, "unique"] & is.na(matches[n, "chr"])) {
                        matches[n, "chr"] = seqname
                        matches[n, "start"] = start(sindex[[n]])
                        matches[n, "strand"] = s
                    } else {
                        matches[n,] = c(NA, NA, NA, FALSE)
                    }
                } else if (length(sindex[[n]]) > 1) {
                    matches[n,] = c(NA, NA, NA, FALSE)
                }
            }
        }
    }

    if (ensemblNotation) {
        matches$chr = gsub("chr", "", matches$chr)
    }

    return(AlignedRead(object,
        id=BStringSet(names(object)),
        chromosome=matches$chr,
        position=matches$start,
        strand=factor(matches$strand)))
}


.alignShortReads_AVASet <- function(object, bsGenome, seqNames, ensemblNotation) {

   aln = alignShortReads(sread(referenceSequences(object)),
       bsGenome, seqNames, ensemblNotation)

   if (!all(names(sread(referenceSequences(object))) == names(sread(aln)))) {
       stop("Some sequences do not match.")
   }
   
   referenceSequences(object)@position = position(aln)
   referenceSequences(object)@chromosome = chromosome(aln)
   referenceSequences(object)@strand = strand(aln)
   
   return(object)
}



setMethod("alignShortReads", 
    signature=signature(object="DNAStringSet", bsGenome="BSgenome",
        seqNames="character", ensemblNotation="logical"),
    .alignShortReads_DNAStringSet)

setMethod("alignShortReads", 
    signature=signature(object="AVASet", bsGenome="BSgenome",
        seqNames="character", ensemblNotation="logical"),
    .alignShortReads_AVASet)

setMethod("alignShortReads", 
    signature=signature(object="DNAStringSet", bsGenome="BSgenome",
        seqNames="character", ensemblNotation="missing"),
    function(object, bsGenome, seqNames) {
        return(.alignShortReads_DNAStringSet(object, bsGenome, seqNames, FALSE))
    })

setMethod("alignShortReads", 
    signature=signature(object="AVASet", bsGenome="BSgenome",
        seqNames="character", ensemblNotation="missing"),
    function(object, bsGenome, seqNames) {
        return(.alignShortReads_AVASet(object, bsGenome, seqNames, FALSE))
    })

setMethod("alignShortReads", 
    signature=signature(object="DNAStringSet", bsGenome="BSgenome",
        seqNames="missing", ensemblNotation="logical"),
    function(object, bsGenome, ensemblNotation) {
        return(.alignShortReads_DNAStringSet(object, bsGenome, names(bsGenome), ensemblNotation))
    })

setMethod("alignShortReads", 
    signature=signature(object="AVASet", bsGenome="BSgenome",
        seqNames="missing", ensemblNotation="logical"),
    function(object, bsGenome, ensemblNotation) {
        return(.alignShortReads_AVASet(object, bsGenome, names(bsGenome), ensemblNotation))
    })

setMethod("alignShortReads", 
    signature=signature(object="DNAStringSet", bsGenome="BSgenome",
        seqNames="missing", ensemblNotation="missing"),
    function(object, bsGenome) {
        return(.alignShortReads_DNAStringSet(object, bsGenome, names(bsGenome), FALSE))
    })

setMethod("alignShortReads", 
    signature=signature(object="AVASet", bsGenome="BSgenome",
        seqNames="missing", ensemblNotation="missing"),
    function(object, bsGenome) {
        return(.alignShortReads_AVASet(object, bsGenome, names(bsGenome), FALSE))
    })






#getAmpliconPositions <- function(dnaSet, genes, dataset="hsapiens_gene_ensembl") {
#  
#  ensembl = useMart("ensembl", dataset=dataset)
#  geneInfo = getBM(attributes=c("hgnc_symbol", "ensembl_gene_id",
#                   "start_position", "end_position", "strand", "chromosome_name"),
#                   filter="hgnc_symbol", values=unique(genes), mart=ensembl)
#
#  geneSeqs = getSequence(id=geneInfo$ensembl_gene_id, type="ensembl_gene_id",
#                         seqType="gene_exon_intron", mart=ensembl)
#
#  gSeqs = DNAStringSet(geneSeqs$gene_exon_intron)
#  names(gSeqs) = geneSeqs$ensembl_gene_id
#  
#  ampPos = data.frame(amplicon=names(dnaSet), start=NA, end=NA, chr=NA,
#              stringsAsFactors=FALSE)
#  for (i in 1:nrow(ampPos)) {
#    amp = dnaSet[[ampPos[i, "amplicon"]]]
#    ind = geneInfo$hgnc_symbol == genes[i]
#    if (geneInfo$strand[ind] == -1) {
#      gene = reverseComplement(gSeqs[[geneInfo[ind, "ensembl_gene_id"]]])
#    } else {
#      gene = gSeqs[[geneInfo[ind, "ensembl_gene_id"]]]
#    }
#
#    m = matchPattern(amp, gene)
#    if (length(m) != 1) {
#      stop(paste("Could not match amplicon:", names(dnaSet)[i]))
#    }
#    
#    ampPos[i, "start"] = geneInfo$start_position[ind] + start(m) - 1
#    ampPos[i, "end"] = geneInfo$start_position[ind] + end(m) - 1
#    ampPos[i, "chr"] = geneInfo$chromosome_name[ind]
#    ampPos[i, "symbol"] = geneInfo$hgnc_symbol[ind]
#    ampPos[i, "ensemblid"] = geneInfo$ensembl_gene_id[ind]
#  }
#
#  return(list(ampPos=ampPos, gSeqs=gSeqs))
#
#}
