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

setGeneric("read454Metrics", function(rochePath, ...) {
    standardGeneric("read454Metrics")})


setGeneric("annotateVariants", function(object, bsGenome) {
    standardGeneric("annotateVariants")})


## generics for plotting functions

setGeneric("plotVariants", signature=c("varData","transcript"),
    function(varData, transcript, legend=TRUE, regions=c(), regLabel="(region label)", col=c("grey","white","black","white"), title="")
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

