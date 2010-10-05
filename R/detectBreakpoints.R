.detectBreakpoints <- function(chimericReads, bpDist=100, minClusterSize=4, removeSoftClips=TRUE, bsGenome) {

    # use library BSgenome.Hsapiens.UCSC.hg19 as default genome for common breakpoint detection
    if(missing(bsGenome)){
      if(!require(BSgenome.Hsapiens.UCSC.hg19))
        stop("library BSgenome.Hsapiens.UCSC.hg19 not found (please try to pass your own genome via function parameters)")
      bsGenome = Hsapiens
    }
    
    # check sam/bam input
    reqFields = c("qname", "rname", "pos", "cigar", "seq", "strand", "qwidth")
    chimericReads[["log"]] = NULL # remove optional filtering log entry
    reads = lapply(names(chimericReads[[1]]), function(elt) {
        if (is.factor(chimericReads[[1]][[elt]])) {
            factor(do.call(c, unname(lapply(chimericReads, function (x) {
                as.character(x[[elt]])
            }))))
                                            
        } else
            do.call(c, unname(lapply(chimericReads, "[[", elt)))
        }
    )
    names(reads) = names(chimericReads[[1]])
    if (!all(is.element(reqFields, names(reads)))) {
        stop(paste("The following SAM fields are necessary:", reqFields))
    }
    reads = do.call("DataFrame", reads)
    reads = reads[order(reads$qname), ]
    
    # remove soft clips, which may be found at the beginning/end of some reads
    if(removeSoftClips){
      message("Removing soft-clipped bases ... ", appendLF = FALSE)
      reads = .removeSoftclips(reads)
      message("done")
    }

    # Build a data frame with one line for each chimeric reads. We can
    # assume that we only have chimeric reads with exactly two local
    # alignments. The data frame has these columns:
    #   p5Chr - Chromosome that the 5' end aligns to.
    #   p5Start - Start position of the 5' alignment. (ref. coordinates)
    #   p5End - End position of the 5' alignment. (reference coordinates)
    #   p5Strand - Strand of the 5' alignment.
    #   p5LocalStart - Start position of the 5' alignment. (read coordinates)
    #   p5LocalEnd - End position of the 5' alignment. (read coordinates)
    #   p5Breakpoint - Breakpoint of the 5' alignment
    #   p3Chr, p3Start, p3End, p3LocalStart, p3LocalEnd, p3Breakpoint - same for 3'
    #
    # Notes: For local coordinates, the first base of the 5' end has
    # coordinate 1. In SAM files, sequences that were aligned to the minus
    # strand of the reference are stored as reverse complement. In that case,
    # the local coordinates do not directly apply to the sequence in the sam
    # file but to its reverse complement.

    crList = split(reads, reads$qname)
    if (!all(sapply(crList, nrow) == 2)) {
        stop("Some of the given reads do not have exactly two
        local alignments.")
    }
    cReads = lapply(crList, function(x) {

        # extract number of matches, dels and ins from cigar
        cigars = extendedCIGARToList(x$cigar)
        dels = sapply(cigars, function(x) sum(x[names(x)=="D"]))
        ins = sapply(cigars, function(x) sum(x[names(x)=="I"]))
        matches = sapply(cigars, function(x) sum(x[names(x)=="M"]))

        # compute global start and end positions
        chr = as.character(x$rname)
        start = x$pos
        end = as.integer(start + matches + dels - 1)

        # compute local start and end position
        # counting starts at the reads 5 prime end
        localStart = as.integer(sapply(cigars, function(cgr) {
            ind = min(which(names(cgr) == "M"))
            cgr = c(0, cgr)
            return(sum(cgr[1:ind]) + 1)
        }))
        localEnd = as.integer(sapply(cigars, function(cgr) {
            ind = max(which(names(cgr) == "M"))
            return(sum(cgr[1:ind]))
        }) - dels)
        ind = x$strand == "-"
        if (sum(ind) > 0) {
            tmp = localStart
            localStart[ind] = x$qwidth[ind] - localEnd[ind] + 1
            localEnd[ind] = x$qwidth[ind] - tmp[ind] + 1
            x$seq[ind] = reverseComplement(x$seq[ind])
            rm(tmp)
        }

        # the smallest local start position is the 5' alignment
        ind5prime = localStart == min(localStart)

        return(data.frame(
            read=x$qname[1],
            p5Chr=chr[ind5prime],
            p5Start=start[ind5prime],
            p5End=end[ind5prime],
            p5Strand=as.character(x$strand[ind5prime]),
            p5LocalStart=localStart[ind5prime],
            p5LocalEnd=localEnd[ind5prime],
            p3Chr=chr[!ind5prime],
            p3Start=start[!ind5prime],
            p3End=end[!ind5prime],
            p3Strand=as.character(x$strand[!ind5prime]),
            p3LocalStart=localStart[!ind5prime],
            p3LocalEnd=localEnd[!ind5prime],
            seq=toString(x$seq[ind5prime]),
            stringsAsFactors=FALSE))
    })
    crInfo = data.frame()
    for (i in 1:length(cReads)) {
        crInfo = rbind(crInfo, cReads[[i]])
    }
    rm(cReads)
    rownames(crInfo) = crInfo$read

    # add breakpoints
    crInfo$p5Breakpoint = NA
    crInfo$p3Breakpoint = NA
    # 5 prime: positive strand
    ind = crInfo$p5Strand == "+"
    crInfo$p5Breakpoint[ind] = crInfo$p5End[ind]
    # 5 prime: negative strand
    ind = crInfo$p5Strand == "-"
    crInfo$p5Breakpoint[ind] = crInfo$p5Start[ind]
    # 3 prime: positive strand
    ind = crInfo$p3Strand == "+"
    crInfo$p3Breakpoint[ind] = crInfo$p3Start[ind]
    # 3 prime: negative strand
    ind = crInfo$p3Strand == "-"
    crInfo$p3Breakpoint[ind] = crInfo$p3End[ind]


    
    # cluster chimeric reads
    dMat = matrix(NA, nrow=nrow(crInfo), ncol=nrow(crInfo))
    rownames(dMat) = rownames(crInfo)
    colnames(dMat) = rownames(crInfo)
    for (i in 1:nrow(dMat)) {
        for (j in i:ncol(dMat)) {

            # first, consider strands: 
            # case 1: both alignment of both reads are from the same strand
            if (crInfo$p5Strand[i] == crInfo$p3Strand[i] &
                crInfo$p5Strand[j] == crInfo$p3Strand[j]) {

                # case 1.1: if the 5' alignments of both reads have the
                #   same strand, both must have the same chromosome
                #   -> same case, same strand
                if (crInfo$p5Strand[i] == crInfo$p5Strand[j] &
                    crInfo$p3Strand[i] == crInfo$p3Strand[j]) {

                    if (crInfo$p5Chr[i] == crInfo$p5Chr[j] &
                        crInfo$p3Chr[i] == crInfo$p3Chr[j]) {

                        dMat[i,j] = (crInfo$p5Breakpoint[i] - crInfo$p5Breakpoint[j])^2 +
                            (crInfo$p3Breakpoint[i] - crInfo$p3Breakpoint[j])^2

                    } else {
                        dMat[i,j] = Inf
                    }

                # case 1.2: if the 5' alignments of both reads are from
                #   different strand, the chromosome of the 5' alignment
                #   of one read must equal the chromosome of the 3'
                #   alignment of the other read
                #   -> same case, different strands
                } else if (crInfo$p5Strand[i] != crInfo$p5Strand[j] &
                    crInfo$p3Strand[i] != crInfo$p3Strand[j]) {

                    if (crInfo$p5Chr[i] == crInfo$p3Chr[j] &
                        crInfo$p3Chr[i] == crInfo$p5Chr[j]) {

                        dMat[i,j] = (crInfo$p5Breakpoint[i] - crInfo$p3Breakpoint[j])^2 +
                            (crInfo$p3Breakpoint[i] - crInfo$p5Breakpoint[j])^2

                    } else {
                        dMat[i,j] = Inf
                    }
                      

                }

            # case 2: alignments of both reads are from different strands
            } else if (crInfo$p5Strand[i] != crInfo$p3Strand[i] &
                crInfo$p5Strand[j] != crInfo$p3Strand[j]) {

                # In this case (strand change), both reads must have the
                # strands for the 5' and 3' alignment.
                if (crInfo$p5Strand[i] == crInfo$p5Strand[j] &
                    crInfo$p3Strand[i] == crInfo$p3Strand[j]) {

                    # case 2.0: Inversions. We know that both reads
                    #   are from the same case, but we cannot
                    #   distinguish between the strands
                    if (crInfo$p5Chr[i] == crInfo$p3Chr[i] &
                        crInfo$p5Chr[j] == crInfo$p3Chr[j] &
                        crInfo$p5Chr[i] == crInfo$p5Chr[j]) {

                        # same strand
                        d1 = (crInfo$p5Breakpoint[i] - crInfo$p5Breakpoint[j])^2 +
                            (crInfo$p3Breakpoint[i] - crInfo$p3Breakpoint[j])^2
                        # different strands
                        d2 = (crInfo$p5Breakpoint[i] - crInfo$p3Breakpoint[j])^2 +
                            (crInfo$p3Breakpoint[i] - crInfo$p5Breakpoint[j])^2
                
                        dMat[i,j] = min(d1, d2)

                    # case 2.1: Translocations. Alignments of both
                    #  reads are from the same chromomsomes and strands
                    #  -> same case, same strand
                    } else if (crInfo$p5Chr[i] == crInfo$p5Chr[j] &
                        crInfo$p3Chr[i] == crInfo$p3Chr[j]) {
                        
                        dMat[i,j] = (crInfo$p5Breakpoint[i] - crInfo$p5Breakpoint[j])^2 +
                            (crInfo$p3Breakpoint[i] - crInfo$p3Breakpoint[j])^2

                    # case 2.2: Translocations. Alignments of both
                    #  reads are from the same strands, but the chr
                    #  of the 5' alignment of one reads equals the
                    #  chr of the 3' alignmnet of the other read.
                    #  -> same case, different strands
                    } else if (crInfo$p5Chr[i] == crInfo$p3Chr[j] &
                        crInfo$p3Chr[i] == crInfo$p5Chr[j]) {

                        dMat[i,j] = (crInfo$p5Breakpoint[i] - crInfo$p3Breakpoint[j])^2 +
                            (crInfo$p3Breakpoint[i] - crInfo$p5Breakpoint[j])^2

                    # different chromosomes
                    } else {
                        dMat[i,j] = Inf
                    }

                
                # Reads are not from the same case
                } else {
                    dMat[i,j] = Inf
                }

            # case 3: both alignments of one read are from the same strand,
            #   whereas the aligment of the other reads are from different
            #   strands -> reads originate from different variants
            } else {
                dMat[i,j] = Inf
            }
        }
    }
    dMat[lower.tri(dMat, diag=FALSE)] = t(dMat)[lower.tri(dMat, diag=FALSE)]
    dMat[dMat == Inf] = max(10^10, 2 * bpDist)

    if (ncol(dMat) > 1) {
        hc = hclust(as.dist(dMat), method="complete")
        bpCluster = cutree(hc, h=bpDist)
        bpClusterSize = sort(table(bpCluster), decreasing=TRUE)
    } else {
        bpCluster = 1
        names(bpCluster) = colnames(dMat)
        bpClusterSize = 1
        names(bpClusterSize) = "1"
    }

    ## only use those cluster for further analysis whose size is greater or equal "minClusterSize" (see function arguments)
    bpClusterSize = bpClusterSize[bpClusterSize >= minClusterSize]
    cluster = names(bpClusterSize)

    ## take care, that reads in "crInfo" and "reads" are in the same order
    crInfo = crInfo[unique(reads$qname), ]
    
    message("Determining consensus breakpoints ... ", appendLF = FALSE)

    ## tolerance: bases read beyond the start and end of the reference sequence
    ## minimum tolerance for single read cluster
    baseTol = 250
    minTol = 10
  
    ## commonBps stores later all information about common Breakpoints and the alignments to the common reference sequences
    seqsC1 = list()
    seqsC2 = list()
    commonBpsC1 = list()
    commonBpsC2 = list()
    commonAlignC1 = list()
    commonAlignC2 = list()
    alignedReadsC1 = list()
    alignedReadsC2 = list()
    
    ## for each cluster, determine a breakpoint, that is shared by all reads (sort cluster from large to small)
    for(clusterInd in cluster){

      clReads = names(bpCluster[bpCluster == clusterInd])
      alnData = crInfo[clReads, ]
      clReads = reads[reads$qname %in% clReads, ]

      numReads = nrow(alnData)
      refSeqs = vector(mode="character")
      localBreakpoints = list()
      brps = list()
      chrs = list()
      strands = list()

      ## determine all possible breakpoints for every read (more than one in case of overlapping read parts)
      start = vector(length=(numReads*2), mode="numeric")
      end = vector(length=(numReads*2), mode="numeric")
      lStart = vector(length=(numReads*2), mode="numeric")
      lEnd = vector(length=(numReads*2), mode="numeric")
      breakpoint = vector(length=(numReads*2), mode="numeric")      

      ## use less tolerance left and right of the reference sequence for single read cluster
      if(numReads > 1)
        tol = baseTol
      else
        tol = minTol

      for(i in 1:numReads){

        
        ## read local, global and breakpoint coordinates (add it to the SAM information in "reads" below)
        if(alnData$p5Start[i] == clReads$pos[2*i-1])
          ind = c(2*i-1, 2*i)
        else
          ind = c(2*i, 2*i-1)          
        start[ind] = c(alnData$p5Start[i], alnData$p3Start[i])
        end[ind] = c(alnData$p5End[i], alnData$p3End[i])
        lStart[ind] = c(alnData$p5LocalStart[i], alnData$p3LocalStart[i])
        lEnd[ind] = c(alnData$p5LocalEnd[i], alnData$p3LocalEnd[i])
        breakpoint[ind] = c(alnData$p5Breakpoint[i], alnData$p3Breakpoint[i])

        ## see if both parts are overlapping or gapped
        seqDiff = alnData$p5LocalEnd[i] - alnData$p3LocalStart[i] + 1

        ## sequences overlap
        if(seqDiff >= 0){
          breakpoints1 = 0:seqDiff
          breakpoints2 = seqDiff:0          
        ## no overlapping (maybe gap between sequences)
        }else{
          breakpoints1 = 0
          breakpoints2 = 0	   	
        }
        
        strand1 = as.character(alnData$p5Strand[i])
        strand2 = as.character(alnData$p3Strand[i])
        chr1 = alnData$p5Chr[i]
        chr2 = alnData$p3Chr[i]
        
        bsGenomeList = list(bsGenome[[paste("chr", chr1, sep="")]], bsGenome[[paste("chr", chr2, sep="")]])
        
        for(j in 1:length(breakpoints1)){

          ## first part of the reference sequence (reverse sequence for neg. strands)
          if(strand1 == "-"){
            refStart1 = alnData$p5Start[i] + breakpoints1[j]
            refEnd1 = alnData$p5End[i] + tol
            refSeqPart1 = toString(subseq(bsGenomeList[[1]], refStart1, refEnd1))
            refSeqPart1 = toString(reverseComplement(DNAString(refSeqPart1)))
            brp1 = refStart1 
          }else{
            refStart1 = alnData$p5Start[i] - tol
            refEnd1 = alnData$p5End[i] - breakpoints1[j]
            refSeqPart1 = toString(subseq(bsGenomeList[[1]], refStart1, refEnd1))
            brp1 = refEnd1
          }
          
          ## second part of the reference sequence (reverse sequence for neg. strands)
          if(strand2 == "-"){
            refStart2 = alnData$p3Start[i] - tol
            refEnd2 = alnData$p3End[i] - breakpoints2[j]
            refSeqPart2 = toString(subseq(bsGenomeList[[2]], refStart2, refEnd2))
            refSeqPart2 = toString(reverseComplement(DNAString(refSeqPart2)))
            brp2 = refEnd2
          }else{
            refStart2 = alnData$p3Start[i] + breakpoints2[j]
            refEnd2 = alnData$p3End[i] + tol
            refSeqPart2 = toString(subseq(bsGenomeList[[2]], refStart2, refEnd2))
            brp2 = refStart2
          }
          
          ## put reference sequence together (first part, gap, second part)
          ## and remember its global and local coordinates, strands, chromosomes
          gap = substr(alnData$seq[i], alnData$p5LocalEnd[i] + 1,  alnData$p3LocalStart[i] - 1)
          brps[[i]] = c(brp1, brp2)
          localBreakpoints[[i]] = c(nchar(refSeqPart1), nchar(refSeqPart1) + nchar(gap) + 1)
          chrs[[i]] = c(chr1, chr2)
          strands[[i]] = c(strand1, strand2)
          refSeqs[i] = paste(refSeqPart1, gap, refSeqPart2, sep="")
        }
      }
      
      ## align reads to reference sequences
      numRefSeqs = length(refSeqs)
      scores = vector(length=numRefSeqs)
      alignments = list(length=numRefSeqs)
      alignedReads = list(length=numRefSeqs)

      for(r in 1:numRefSeqs){
        
        refSeq = refSeqs[r]
        refBps = brps[[r]]
        refp5Chr = chrs[[r]]
        alnReads = vector(length=numReads)
        alnChrs = vector(length(numReads), mode="character")
        alnStrands = vector(length(numReads), mode="character")
        
        ## use the reverse complement of sequences on complement strand
        ## (i.e. 5' chromosome differs between read and reference,
        ## or (in case of inversions) p5Breakpoint differs from 5' reference breakpoint:
        ## distance between breakpoints > 1.000bp (empirical threshold))
        for(j in 1:numReads){
          if(alnData$p5Chr[j] != refp5Chr[1] | abs(alnData$p5Breakpoint[j] - refBps[1]) > 1000){
            alnChrs[j] = paste(alnData$p3Chr[j], alnData$p5Chr[j],sep="/")
            alnStrands[j] = paste(alnData$p3Strand[j], alnData$p5Strand[j],sep="/")
            alnReads[j] = toString(reverseComplement(DNAString(alnData$seq[j])))
          }else{
            alnChrs[j] = paste(alnData$p5Chr[j], alnData$p3Chr[j],sep="/")
            alnStrands[j] = paste(alnData$p5Strand[j], alnData$p3Strand[j],sep="/")
            alnReads[j] = toString(alnData$seq[j])
          }
        }
        alnReads = DNAStringSet(alnReads)
        names(alnReads) = alnData$read
        
        ## alignment
        substitutionMatrix = nucleotideSubstitutionMatrix(match=1, mismatch=-1)
        alignment = pairwiseAlignment(pattern=alnReads, subject=refSeq, 
          type="global-local", substitutionMatrix=substitutionMatrix, gapOpening=0, gapExtension=-1)
        scores[r] = sum(score(alignment))
        alignments[[r]] = alignment
        alignedReads[[r]] = AlignedRead(sread=alnReads, chromosome=factor(alnChrs), position=start(Views(alignment)), strand=factor(alnStrands))
      
      }
      ## add local, global and breakpoint coordinates to the alignment information in "reads"
      clReads$start = as.integer(start)
      clReads$end = as.integer(end)
      clReads$localStart = as.integer(lStart)
      clReads$localEnd = as.integer(lEnd)
      clReads$breakpoint = as.integer(breakpoint)
      rm(list=c("start", "end", "lStart", "lEnd", "breakpoint"))
      
      ## select breakpoints with maximum score, and the corresponding alignment
      maxScoreInd = which(scores == max(scores))[1]
      seqsC1[[clusterInd]] = alnData
      seqsC2[[clusterInd]] = data.frame()
      commonBpsC1[[clusterInd]] = data.frame(
                   "chr" = chrs[[maxScoreInd]],
                   "strand" = strands[[maxScoreInd]],
                   "breakpoint" = brps[[maxScoreInd]],
                   "localStart" = c(1, localBreakpoints[[maxScoreInd]][2]),
                   "localEnd" = c(localBreakpoints[[maxScoreInd]][1], nchar(refSeqs[maxScoreInd])),
                   "refSeq" = refSeqs[maxScoreInd],
                   stringsAsFactors = FALSE)
      commonBpsC2[[clusterInd]] = data.frame()
      commonAlignC1[[clusterInd]] = alignments[[maxScoreInd]]
      commonAlignC2[[clusterInd]] = NA
      alignedReadsC1[[clusterInd]] = alignedReads[[maxScoreInd]]
      alignedReadsC2[[clusterInd]] = AlignedRead()

    }
    if(length(cluster >= 1)){
      bpNames = paste("BP", 1:length(cluster), sep="")
      names(seqsC1) = bpNames
      names(seqsC2) = bpNames
      names(commonBpsC1) = bpNames
      names(commonBpsC2) = bpNames
      names(commonAlignC1) = bpNames
      names(commonAlignC2) = bpNames
      names(alignedReadsC1) = bpNames
      names(alignedReadsC2) = bpNames
    }    

    commonBps = new("Breakpoints", seqsC1=seqsC1, commonBpsC1=commonBpsC1, commonAlignC1=commonAlignC1, alignedReadsC1=alignedReadsC1, seqsC2=seqsC2, commonBpsC2=commonBpsC2, commonAlignC2=commonAlignC2, alignedReadsC2=alignedReadsC2)
    message("done")
    
    message("Detected ", length(cluster), " cluster of size >= ", minClusterSize)
    
    return(commonBps)

    
}

setMethod("detectBreakpoints",
          signature=signature(chimericReads="list"),
          .detectBreakpoints)


## Methods for preprocessing:


## Some reads may contain soft-clipped bases at the beginning of the first part of the read or 
## at the end of the second part. This function removes these unaligned (sub-)sequences
## and adjusts the cigar string, the sequence, the sequence width (qwidth) and the 
## local start/end coordinates

.removeSoftclips <- function(seqs){

   ## initialize vectors for updated entries (for some reason, updating entries directly in the IRanges DataFrame is relatively slow)
   numSeqs = nrow(seqs)
   newCigar = vector(length=numSeqs, mode="character")
   newSeq = vector(length=numSeqs, mode="character")
   newWidth = vector(length=numSeqs, mode="numeric")
   ## keep sequence names (if existent)
   if(!is.null(names(seqs$seq)))
      names(newSeq) = names(seqs$seq)

   ## 1. compute local start from cigar string and convert it for neg. strands

   # extract number of matches, dels and ins from cigar
   cigars = extendedCIGARToList(seqs$cigar)
   dels = sapply(cigars, function(x) sum(x[names(x)=="D"]))
   ins = sapply(cigars, function(x) sum(x[names(x)=="I"]))
   matches = sapply(cigars, function(x) sum(x[names(x)=="M"]))

   # compute local start and end position
   # counting starts at the reads 5 prime end
   localStart = as.integer(sapply(cigars, function(cgr) {
       ind = min(which(names(cgr) == "M"))
       cgr = c(0, cgr)
       return(sum(cgr[1:ind]) + 1)
   }))
   localEnd = as.integer(sapply(cigars, function(cgr) {
       ind = max(which(names(cgr) == "M"))
       return(sum(cgr[1:ind]))
   }) - dels)
   ind = seqs$strand == "-"
   if (sum(ind) > 0) {
       tmp = localStart
       localStart[ind] = seqs$qwidth[ind] - localEnd[ind] + 1
       localEnd[ind] = seqs$qwidth[ind] - tmp[ind] + 1
       rm(tmp)
   }
   
   for(i in seq(1, numSeqs, 2)){

     lStart1 = localStart[i]
     lStart2 = localStart[i + 1]

     ## 2. find out order of reads
     if(lStart1 < lStart2){
       part1 = i
       part2 = i + 1
     }else{
       part1 = i + 1
       part2 = i
     }

     ## 3. remove soft clip from the cigar string and the sequence itself
     ##    unify reads by reversing cigar string and sequence for reads on the negative strand
     ##    -> soft clip is the the beginning of the first part and at the end of the second part
     cig1 = extendedCIGARToList(seqs$cigar[part1])[[1]]
     cig2 = extendedCIGARToList(seqs$cigar[part2])[[1]]
     seq1 = seqs$seq[part1]
     seq2 = seqs$seq[part2]             
     sc1 = 0
     sc2 = 0

     if(seqs$strand[part1] == "-")
       cig1 = rev(cig1)
     if(seqs$strand[part2] == "-")
       cig2 = rev(cig2)

     ## test if any soft clips exist and remove their entry
     if(names(cig1[1]) == "S"){
       sc1 = cig1[1]               
       cig1 = cig1[-1]
     }
     if(names(cig2[length(cig2)]) == "S"){
       sc2 = cig2[length(cig2)]
       cig2 = cig2[-(length(cig2))]
     }

     ## adjust coordinates and cigar entries if soft clip is existent
     if(sc1 > 0 || sc2 > 0){

       ## adjust other cigar entries before/after the soft clip
       cig1[length(cig1)] = cig1[length(cig1)] - sc2
       cig2[1] = cig2[1] - sc1

       ## remove soft clip from sequence
       if(seqs$strand[part1] == "-")
         seq1 = reverseComplement(seq1)
       if(seqs$strand[part2] == "-")    
         seq2 = reverseComplement(seq2)                              
       seq1 = subseq(seq1, sc1 + 1, width(seq1) - sc2)
       seq2 = seq1
       if(seqs$strand[part1] == "-"){
         cig1 = rev(cig1)
         seq1 = reverseComplement(seq1)
       }
       if(seqs$strand[part2] == "-"){
         cig2 = rev(cig2)                 
         seq2 = reverseComplement(seq2)
       }
       
       ## update breakpoints object
       newCigar[c(part1, part2)] = c(listToExtendedCIGAR(list(cig1))[1], listToExtendedCIGAR(list(cig2))[1]) 
       newSeq[c(part1, part2)] = c(as.character(seq1), as.character(seq2))
       newWidth[c(part1, part2)] = c(seqs$qwidth[part1] - sc1 - sc2, seqs$qwidth[part2] - sc1 - sc2)
       
     }else{
       ## if no soft clip exists, keep old coordinates and cigar strings
       newCigar[c(part1, part2)] = c(seqs$cigar[part1], seqs$cigar[part2])
       newSeq[c(part1, part2)] = c(toString(seqs$seq[part1]), toString(seqs$seq[part2]))
       newWidth[c(part1, part2)] = c(seqs$qwidth[part1], seqs$qwidth[part2])
     }
   }

   ## update given IRanges DataFrame
   seqs$cigar = newCigar
   seqs$seq = DNAStringSet(newSeq, use.names=TRUE)
   seqs$qwidth = as.integer(newWidth)
   return(seqs)
}



## computes local coordinates (start/end of mapped region within a read) from a given cigar 
## note: cigar must be given as vector (see "extendedCIGARToList")

.computeLocalCoordinates <- function(cigar){
  
  dels = sum(cigar[names(cigar)=="D"])

  ## local start
  ind = min(which(names(cigar) == "M"))
  tmp = c(0, cigar)
  localStart = sum(tmp[1:ind]) + 1

  ## local end
  dels = sum(cigar[names(cigar)=="D"])
  ind = max(which(names(cigar) == "M"))
  localEnd = sum(cigar[1:ind]) - dels

  return(c(localStart, localEnd))
}


