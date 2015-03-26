.ava2vcf <- function(object, filename, annot) {
      
  # get the index of the reference for each variant
  refSeqInd <- match(fData(object)$referenceSeqID, as.character(id(referenceSequences(object))))
  
  # get the chromosome and strand for each variant
  chromosome <- as.character(chromosome(referenceSequences(object)[refSeqInd]))
  posStrand <- strand(referenceSequences(object))[refSeqInd] == "+"
  if(sum(is.na(chromosome)) > 0 | sum(is.na(posStrand)) > 0) {
   stop("Chromosome and/or strand information for references is missing, maybe you have to run alignShortReads first.") 
  }

  ### calculate the absolute position on the refrence genome for each variant ###
  
  # calculate the absolute positions of the reference sequences
  startRefSeqs <- position(referenceSequences(object)[refSeqInd])
  endRefSeqs <-  startRefSeqs + width(sread(referenceSequences(object)[refSeqInd])) - 1

  seq_start <- fData(object)$start
  seq_end <- fData(object)$end
  
  abs_start <- numeric(length(chromosome))
  # as insertions have .5 positions we have to correct them
  abs_start[posStrand] <- startRefSeqs[posStrand] + ceiling(seq_start[posStrand]) - 1  
  abs_start[!posStrand] <- endRefSeqs[!posStrand] - floor(seq_end[!posStrand]) + 1
  
  # not needed as vcf only uses start positions
  # abs_end = numeric(length(chromosome))
  # as insertions have .5 positions we have to correct them
  # abs_end[posStrand] <- startRefSeqs[posStrand] + ceiling(seq_end[posStrand]) - 1
  # abs_end[!posStrand] <- endRefSeqs[!posStrand] - floor(seq_start[!posStrand]) + 1
  
  mutation_type <- sapply(fData(object)$canonicalPattern, function(x) {substr(x, 1, 1)})
  mutation_type[mutation_type == "s"] <- "point"
  mutation_type[mutation_type == "d"] <- "deletion"
  mutation_type[mutation_type == "i"] <- "insertion"
  
  #get reference and alt alleles
  ref <- fData(object)$referenceBases
  alt <- fData(object)$variantBase
  
  # plus strand
  
  # nothing to do for point mutations on plus strand

  # add base before deletion on plus strand for deletions
  ind <- mutation_type == "deletion" & posStrand
  if(sum(ind) > 0) {
    ref[ind] <- paste(getSeq(Hsapiens, names=chromosome[ind], start=abs_start[ind]-1, width=1, as.character=TRUE), ref[ind], sep="")
    alt[ind] <- getSeq(Hsapiens, names=chromosome[ind], start=abs_start[ind]-1, width=1, as.character=TRUE)
    # new start position is that of preceding base
    abs_start[ind] <- abs_start[ind] - 1
  }

  # add base before insertion on plus strand for insertions
  ind <- mutation_type == "insertion" & posStrand
  if(sum(ind) > 0) {
    ref[ind] <- getSeq(Hsapiens, names=chromosome[ind], start=abs_start[ind]-1, width=1, as.character=TRUE)
    alt[ind] <- paste(getSeq(Hsapiens, names=chromosome[ind], start=abs_start[ind]-1, width=1, as.character=TRUE), alt[ind], sep="")
    # new start position is that of preceding base
    abs_start[ind] <- abs_start[ind] - 1
  }
 
  # minus strand
  # for mutations with a reference on the minus strand we have to change bases and directions
  
  # invert point mutations on minus strand
  ind <- mutation_type == "point" & !posStrand
  if(sum(ind) > 0) {
    ref[ind] <- as.character(reverseComplement(DNAStringSet(ref[ind])))
    alt[ind] <- as.character(reverseComplement(DNAStringSet(alt[ind])))
  }

  # add base before deletion and create reverse complement
  ind <- mutation_type == "deletion" & !posStrand
  if(sum(ind) > 0) {
    ref[ind] <- paste(getSeq(Hsapiens, names=chromosome[ind], start=abs_start[ind]-1, width=1, as.character=TRUE), as.character(reverseComplement(DNAStringSet(ref[ind]))), sep="")
    alt[ind] <- getSeq(Hsapiens, names=chromosome[ind], start=abs_start[ind]-1, width=1, as.character=TRUE)
    # new end position is that of preceding base
    abs_start[ind] <- abs_start[ind] - 1
  }

  # add base before (after in + coordinates) insertion on minus strand and reverse complement
  ind <- mutation_type == "insertion" & !posStrand
  if(sum(ind) > 0) {
    ref[ind] <- getSeq(Hsapiens, names=chromosome[ind], start=abs_start[ind]-1, width=1, as.character=TRUE)
    alt[ind] <- paste(getSeq(Hsapiens, names=chromosome[ind], start=abs_start[ind]-1, width=1, as.character=TRUE), as.character(reverseComplement(DNAStringSet(alt[ind]))), sep="")
    # new end position is that of preceding base
    abs_start[ind] <- abs_start[ind] - 1
  }
  
  # get known SNPs
  if(!missing(annot)) {
    anna <- annotatedVariants(annot)
    snps <- sapply(anna, function(x)
        if(ncol(x$snps) > 0) {
          # only recognize variants with exactly one known SNP
          rss <- grep("rs", x$snps$refsnp_id, fixed=TRUE)
          if(length(rss == 1)) {
            as.character(x$snps$refsnp_id[rss])
          } else {
            "."
          }
        } else {
          "."
        }
    )
    snpflags <- !(snps == ".")

  }
  
  ### create vcf object ###
  
  # header
  h_meta <- DataFrame(Value=c("VCFv4.1", gsub("-", "", Sys.Date(), fixed=TRUE), "R453Plus1Toolbox", "BSgenome.Hsapiens.UCSC.hg19"), 
                      row.names=c("fileFormat", "fileDate", "source", "reference"))
  if(!missing(annot)) {
    h_info <- DataFrame(Number=c("1", "4", "A", "A", "0"), Type=c("Integer", "Integer", "Float", "Integer", "Flag"), 
                        Description=c("Total Depth", 
                                      "Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles",
                                      "Allele Frequency, for each ALT allele, in the same order as listed",
                                      "Allele count in genotypes, for each ALT allele, in the same order as listed",
                                      "dbSNP membership, build 135")) # INFO zur variante, welche sind möglich
    rownames(h_info) <- c("DP", "DP4", "AF", "AC", "DB")
  } else {
    h_info <- DataFrame(Number=c("1", "4", "A", "A"), Type=c("Integer", "Integer", "Float", "Integer"), 
                        Description=c("Total Depth", 
                                      "Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles",
                                      "Allele Frequency, for each ALT allele, in the same order as listed",
                                      "Allele count in genotypes, for each ALT allele, in the same order as listed"
                        )) # INFO zur variante, welche sind möglich
    rownames(h_info) <- c("DP", "DP4", "AF", "AC")
  }

  h_format <- DataFrame(Number=c("1", "."), Type=c("Integer", "Integer"), 
                        Description=c("Read Depth", "Allelic depths for the alt alleles in the order listed")) #FORMAT zu den einzelnen samples
  rownames(h_format) <- c("DP", "AD")
  
  header <- VCFHeader(reference="BSgenome.Hsapiens.UCSC.hg19", samples=as.character(pData(object)$SampleID), 
                      header=DataFrameList(META=h_meta, INFO=h_info, FORMAT=h_format))
  exptData <- SimpleList(header=header)

  # body
  rowData <- GRanges(seqnames=Rle(chromosome), ranges=IRanges(start=abs_start, end=abs_start), strand="+", paramRangeID=rownames(fData(object)))
  if(!missing(annot)) {
    names(rowData) <- snps
  } else {
    names(rowData) <- rep(".", length(rowData))
  }
        
  colData <- DataFrame(Samples=pData(object)$Annotation, row.names=pData(object)$SampleID)

  if(!missing(annot)) {
    info <- DataFrame(DP=rowSums(assayData(object)$totalForwCount + assayData(object)$totalRevCount),
                      DP4=paste(rowSums(assayData(object)$totalForwCount - assayData(object)$variantForwCount),
                                rowSums(assayData(object)$totalRevCount - assayData(object)$variantRevCount),
                                rowSums(assayData(object)$variantForwCount),
                                rowSums(assayData(object)$variantRevCount),
                                sep=","),
                      AF=round(rowSums(assayData(object)$variantForwCount + assayData(object)$variantRevCount) / 
                         rowSums(assayData(object)$totalForwCount + assayData(object)$totalRevCount), 4),
                      AC=rowSums(assayData(object)$variantForwCount + assayData(object)$variantRevCount),
                      DB=snpflags,
                      row.names=rownames(fData(object)))
  } else {
    info <- DataFrame(DP=rowSums(assayData(object)$totalForwCount + assayData(object)$totalRevCount),
                      DP4=paste(rowSums(assayData(object)$totalForwCount - assayData(object)$variantForwCount),
                                rowSums(assayData(object)$totalRevCount - assayData(object)$variantRevCount),
                                rowSums(assayData(object)$variantForwCount),
                                rowSums(assayData(object)$variantRevCount),
                                sep=","),
                      AF=round(rowSums(assayData(object)$variantForwCount + assayData(object)$variantRevCount) / 
                         rowSums(assayData(object)$totalForwCount + assayData(object)$totalRevCount), 4),
                      AC=rowSums(assayData(object)$variantForwCount + assayData(object)$variantRevCount),
                      row.names=rownames(fData(object)))
  }
  
  format <- SimpleList(DP=(assayData(object)$totalForwCount + assayData(object)$totalRevCount), 
                       AD=(assayData(object)$variantForwCount + assayData(object)$variantRevCount))
  
  fixed <- DataFrame(REF=DNAStringSet(ref), ALT=CharacterList(as.list(alt)), QUAL=0, FILTER=as.character("."))
  
  # merge it all together
  vcf <- VCF(rowRanges=rowData, colData=colData, exptData=exptData, fixed=fixed, geno=format, info=info)

  # return VCF object or write it to file
  if(missing(filename)) {
    return(vcf)
  } else {
    writeVcf(vcf, filename) 
  }
}

setMethod("ava2vcf",
    signature=signature(object="AVASet"),
    .ava2vcf)
