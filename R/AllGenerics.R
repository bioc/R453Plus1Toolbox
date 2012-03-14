# class Breakpoints
setGeneric("seqsC1", function(object) {
    standardGeneric("seqsC1")})
setGeneric("seqsC2", function(object) {
    standardGeneric("seqsC2")})
setGeneric("commonBpsC1", function(object) {
    standardGeneric("commonBpsC1")})
setGeneric("commonBpsC2", function(object) {
    standardGeneric("commonBpsC2")})
setGeneric("commonAlignC1", function(object) {
    standardGeneric("commonAlignC1")})
setGeneric("commonAlignC2", function(object) {
    standardGeneric("commonAlignC2")})
setGeneric("alignedReadsC1", function(object) {
    standardGeneric("alignedReadsC1")})
setGeneric("alignedReadsC2", function(object) {
    standardGeneric("alignedReadsC2")})

setGeneric("seqsC1<-", function(object, value) {
    standardGeneric("seqsC1<-")})
setGeneric("seqsC2<-", function(object, value) {
    standardGeneric("seqsC2<-")})
setGeneric("commonBpsC1<-", function(object, value) {
    standardGeneric("commonBpsC1<-")})
setGeneric("commonBpsC2<-", function(object, value) {
    standardGeneric("commonBpsC2<-")})
setGeneric("commonAlignC1<-", function(object, value) {
    standardGeneric("commonAlignC1<-")})
setGeneric("commonAlignC2<-", function(object, value) {
    standardGeneric("commonAlignC2<-")})
setGeneric("alignedReadsC1<-", function(object, value) {
    standardGeneric("alignedReadsC1<-")})
setGeneric("alignedReadsC2<-", function(object, value) {
    standardGeneric("alignedReadsC2<-")})



setGeneric("read454QualityFilterMetrics", function(object) {
    standardGeneric("read454QualityFilterMetrics")})

setGeneric("dataTotal", function(object) {
    standardGeneric("dataTotal")})

setGeneric("dataTotal<-", function(object, value) {
    standardGeneric("dataTotal<-")})

setGeneric("dataTCAG", function(object) {
    standardGeneric("dataTCAG")})

setGeneric("dataTCAG<-", function(object, value) {
    standardGeneric("dataTCAG<-")})

setGeneric("dataCATG", function(object) {
    standardGeneric("dataCATG")})

setGeneric("dataCATG<-", function(object, value) {
    standardGeneric("dataCATG<-")})

setGeneric("read454BaseCallerMetrics", function(object) {
    standardGeneric("read454BaseCallerMetrics")})

setGeneric("lengthTCAG", function(object) {
    standardGeneric("lengthTCAG")})

setGeneric("lengthTCAG<-", function(object, value) {
    standardGeneric("lengthTCAG<-")})

setGeneric("scoreTCAG", function(object) {
    standardGeneric("scoreTCAG")})

setGeneric("scoreTCAG<-", function(object, value) {
    standardGeneric("scoreTCAG<-")})

setGeneric("lengthCATG", function(object) {
    standardGeneric("lengthCATG")})

setGeneric("lengthCATG<-", function(object, value) {
    standardGeneric("lengthCATG<-")})

setGeneric("scoreCATG", function(object) {
    standardGeneric("scoreCATG")})

setGeneric("scoreCATG<-", function(object, value) {
    standardGeneric("scoreCATG<-")})

setGeneric("annotatedVariants", function(object) {
    standardGeneric("annotatedVariants")})

setGeneric("annotatedVariants<-", function(object, value) {
    standardGeneric("annotatedVariants<-")})

setGeneric("plotVariationFrequency", function(object, plotRange, ...) {
    standardGeneric("plotVariationFrequency")})

setGeneric("alignShortReads", function(object, bsGenome, seqNames, ensemblNotation) {
    standardGeneric("alignShortReads")})



setGeneric("genomeSequencerMIDs", function(mid) {
    standardGeneric("genomeSequencerMIDs")})

setGeneric("sequenceCaptureLinkers", function(name) {
    standardGeneric("sequenceCaptureLinkers")})

setGeneric("demultiplexReads", function(reads, mids, numMismatches, trim) {
    standardGeneric("demultiplexReads")})

setGeneric("removeLinker", function(reads, linker, removeReadsWithoutLinker, minOverlap, penalty) {
    standardGeneric("removeLinker")})

setGeneric("filterChimericReads", function(alnReads, targetRegion, linkerSeq, minDist, dupReadDist) {
    standardGeneric("filterChimericReads")})

setGeneric("detectBreakpoints", signature=c("chimericReads"),
           function(chimericReads, bpDist=100, minClusterSize=4, removeSoftClips=TRUE, bsGenome){
             standardGeneric("detectBreakpoints")})

setGeneric("mergeBreakpoints", function(breakpoints, maxDist, mergeBPs) {
    standardGeneric("mergeBreakpoints")})


setGeneric("readsOnTarget", function(alnReads, targetRegion) {
    standardGeneric("readsOnTarget")})

setGeneric("coverageOnTarget", function(alnReads, targetRegion) {
    standardGeneric("coverageOnTarget")})

setGeneric("calculateTiTv", function(object) {
    standardGeneric("calculateTiTv")})

setGeneric("read454Metrics", function(rochePath, ...) {
    standardGeneric("read454Metrics")})


setGeneric("annotateVariants", function(object, bsGenome) {
    standardGeneric("annotateVariants")})


## generics for plotting functions

setGeneric("plotVariants", signature=c("data", "gene"),
    function(data, gene, transcript=NA, regions, mutationInfo, horiz=FALSE, cex=1, title="", legend=TRUE)
        standardGeneric("plotVariants"))

setGeneric("plotChimericReads", signature=("brpData"),
    function(brpData, geneSymbols=FALSE, plotMut=TRUE, plotBasePairs=FALSE, maxBasePairs=50, legend=FALSE, title="", 
	col=c("red", "green", "black", "orange"))
        standardGeneric("plotChimericReads"))

setGeneric("plotAmpliconCoverage", signature=c("avaSet", "type", "bothDirections"),
     function(avaSet, type="amplicon", bothDirections=TRUE,
         cex.names=0.8, cex.axis=0.8, las=3,
         col=c(rgb(217/255, 214/255, 209/255), rgb(173/255, 38/255, 36/255)),
         ...) {
         standardGeneric("plotAmpliconCoverage")})

setGeneric("htmlReport", signature=c("object"),
    function(object, annot, blocks=c(), transcripts=c(), sampleCols, minMut=3, dir="HTMLReport", title="Summary")
        standardGeneric("htmlReport"))


## generics for filterVariants

setGeneric("getVariantPercentages", signature=c("object"),
    function(object, direction="both")
        standardGeneric("getVariantPercentages"))

setGeneric("setVariantFilter", signature=c("object"),
    function(object, filter=0)
        standardGeneric("setVariantFilter"))

setGeneric("variantFilterPerc", signature = c("object"),
    function(object)
        standardGeneric("variantFilterPerc"))

setGeneric("variantFilter", signature = c("object"),
     function(object)
        standardGeneric("variantFilter"))


## generics for the AVASet:

setGeneric("AVASet", function(dirname) 
    standardGeneric("AVASet"))

setGeneric("readSampleData", function(dir_projectDef) 
    standardGeneric("readSampleData"))

setGeneric("readVariants", function(dir_variants,samples,refSeqs) 
    standardGeneric("readVariants"))

setGeneric("readReferenceSequences", function(dir_projectDef) 
    standardGeneric("readReferenceSequences"))

setGeneric("readAmplicons", function(dir_projectDef, dir_align, samples) 
    standardGeneric("readAmplicons"))

setGeneric("readAmpliconPositions", signature=c("object", "dnaSet","genes"), 
    function(object, dnaSet, genes, dataset="hsapiens_gene_ensembl") 
	standardGeneric("readAmpliconPositions"))

setGeneric("fDataAmp", function(object) 
    standardGeneric("fDataAmp"))

setGeneric("featureDataAmp", function(object) 
    standardGeneric("featureDataAmp"))

setGeneric("assayDataAmp", function(object) 
    standardGeneric("assayDataAmp"))

setGeneric("referenceSequences", signature = c("object"), 
    function(object) 
	standardGeneric("referenceSequences"))

setGeneric("getAlignedReads", signature = c("object"), 
    function(object, amplicons, dir) 
	standardGeneric("getAlignedReads"))

setGeneric("featureDataAmp<-", signature = c("object","value"), 
    function(object, value) 
	standardGeneric("featureDataAmp<-"))

setGeneric("assayDataAmp<-", signature = c("object","value"),
    function(object, value) 
	standardGeneric("assayDataAmp<-"))

setGeneric("referenceSequences<-", signature = c("object","value"), 
    function(object, value) 
	standardGeneric("referenceSequences<-"))


#generics for the MapperSet:

setGeneric("MapperSet", function(dirs, samplenames)
    standardGeneric("MapperSet"))

setGeneric("getReadStatus", function(object)
    standardGeneric("getReadStatus"))



## generics for the AVASet and MapperSet:

setGeneric("dirs", signature = c("object"), 
    function(object) 
	standardGeneric("dirs"))

setGeneric("variantFilterPerc<-", signature = c("object","value"), 
    function(object, value) 
	standardGeneric("variantFilterPerc<-"))

setGeneric("variantFilter<-", signature = c("object","value"), 
    function(object, value) 
	standardGeneric("variantFilter<-"))
	
	
## generics for SFFRead and SFFContainer

setGeneric("name<-", function(x, value) standardGeneric("name<-"))

setGeneric("flowgramFormat", function(x) standardGeneric("flowgramFormat"))
setGeneric("flowgramFormat<-", function(x, value) standardGeneric("flowgramFormat<-"))

setGeneric("flowChars", function(x) standardGeneric("flowChars"))
setGeneric("flowChars<-", function(x, value) standardGeneric("flowChars<-"))

setGeneric("keySequence", function(x) standardGeneric("keySequence"))
setGeneric("keySequence<-", function(x, value) standardGeneric("keySequence<-"))

setGeneric("clipQualityLeft", function(x) standardGeneric("clipQualityLeft"))
setGeneric("clipQualityLeft<-", function(x, value) standardGeneric("clipQualityLeft<-"))

setGeneric("clipQualityRight", function(x) standardGeneric("clipQualityRight"))
setGeneric("clipQualityRight<-", function(x, value) standardGeneric("clipQualityRight<-"))

setGeneric("clipAdapterLeft", function(x) standardGeneric("clipAdapterLeft"))
setGeneric("clipAdapterLeft<-", function(x, value) standardGeneric("clipAdapterLeft<-"))

setGeneric("clipAdapterRight", function(x) standardGeneric("clipAdapterRight"))
setGeneric("clipAdapterRight<-", function(x, value) standardGeneric("clipAdapterRight<-"))

setGeneric("flowgram", function(x) standardGeneric("flowgram"))
setGeneric("flowgram<-", function(x, value) standardGeneric("flowgram<-"))

setGeneric("flowgrams", function(x) standardGeneric("flowgrams"))
setGeneric("flowgrams<-", function(x, value) standardGeneric("flowgrams<-"))

setGeneric("flowIndexes", function(x) standardGeneric("flowIndexes"))
setGeneric("flowIndexes<-", function(x, value) standardGeneric("flowIndexes<-"))

setGeneric("read", function(x) standardGeneric("read"))
setGeneric("read<-", function(x, value) standardGeneric("read<-"))

setGeneric("quality<-", function(x, value) standardGeneric("quality<-"))

setGeneric("reads", function(x) standardGeneric("reads"))
setGeneric("reads<-", function(x, value) standardGeneric("reads<-"))

setGeneric("addRead", function(x, read) standardGeneric("addRead"))
setGeneric("getRead", function(x, readname) standardGeneric("getRead"))

setGeneric(name="sff2fastq", 
           def=function(x, outdir, fname) standardGeneric("sff2fastq"),
           signature=c("x"))


## generics for quality functions
setGeneric("readLengthStats", function(object) standardGeneric("readLengthStats"))
setGeneric(name="readLengthHist", 
           def=function(object, cutoff=0.99, xlab="Read length", ylab="Number of sequences", 
             col="firebrick1", breaks=100, ...) standardGeneric("readLengthHist"),
           signature=c("object"))

setGeneric("baseQualityStats", function(object) standardGeneric("baseQualityStats"))
setGeneric(name="baseQualityHist", 
           def=function(object, xlab="Quality score", ylab="Number of bases", 
             col="firebrick1", breaks=40, ...) standardGeneric("baseQualityHist"),
           signature=c("object"))
setGeneric(name="sequenceQualityHist", 
           def=function(object, xlab="Mean of quality scores per sequence", 
             ylab="Number of sequences", col="firebrick1", ...) 
             standardGeneric("sequenceQualityHist"),
           signature=c("object"))
setGeneric(name="positionQualityBoxplot", 
           def=function(object, range, binsize=10, 
             xlab=paste("Read position in bp (Bin size: ", binsize, "bp)", sep=""), 
             ylab="Quality score", col="firebrick1", ...) standardGeneric("positionQualityBoxplot"),
           signature=c("object"))

setGeneric("baseFrequency", function(object) standardGeneric("baseFrequency"))
setGeneric(name="nucleotideCharts", 
           def=function(object, range=0.95, linetypes=c(A="l", C="l", G="l", T="l", N="l"), 
             linecols=c(A="black", C="red", G="blue", T="green", N="grey70"), xlab="Position",
             ylab="Frequency", ...) standardGeneric("nucleotideCharts"),
           signature=c("object"))
setGeneric("gcContent", function(object) standardGeneric("gcContent"))
setGeneric(name="gcPerPosition", 
           def=function(object, range=0.95, type="l", col="blue", xlab="Position", ylab="Frequency",
             ...) standardGeneric("gcPerPosition"),
           signature=c("object"))
setGeneric(name="gcContentHist", 
           def=function(object, xlab="GC content", ylab="Number of sequences", col="firebrick1",
             breaks=50, ...) standardGeneric("gcContentHist"),
           signature=c("object"))

setGeneric(name="complexity.dust", 
           def=function(object, xlab="Complexity score (0=high, 100=low)", 
             ylab="Number of sequences", xlim=c(0, 100), col="firebrick1", breaks=100, ...) 
             standardGeneric("complexity.dust"),
           signature=c("object"))
setGeneric(name="complexity.entropy", 
           def=function(object, xlab="Complexity score (0=low, 100=high)", 
             ylab="Number of sequences", xlim=c(0, 100), col="firebrick1", breaks=100, ...) 
             standardGeneric("complexity.entropy"),
           signature=c("object"))

setGeneric(name="dinucleotideOddsRatio", 
           def=function(object, xlab="Under-/over-representation of dinucleotides", 
             col="firebrick1", ...) standardGeneric("dinucleotideOddsRatio"),
           signature=c("object"))

setGeneric(name="flowgramBarplot",
           def=function(x, range=c(0, length(flowgram(x))), xlab="Flow sequence", ylab="Flow intensity",
             col=c(A="black", C="red", G="blue", T="green"), ...)
             standardGeneric("flowgramBarplot"),
           signature=c("x"))

setGeneric(name="homopolymerHist",
           def=function(x, range=c(0, length(flowgram(x))), xlab="Homopolymer length",
             ylab="Number of homopolymers", col=c(A="black", C="red", G="blue", T="green"), ...)
             standardGeneric("homopolymerHist"),
           signature=c("x"))
