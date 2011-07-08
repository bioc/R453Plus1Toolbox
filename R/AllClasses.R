setClass("Roche454BaseCallerMetrics",
    representation=list(
        dataTCAG="data.frame",
        dataCATG="data.frame",
        lengthTCAG="matrix",
        scoreTCAG="matrix",
        lengthCATG="matrix",
        scoreCATG="matrix"),
    prototype=prototype(
        dataTCAG=data.frame(),
        dataCATG=data.frame(),
        lengthTCAG=matrix(),
        scoreTCAG=matrix(),
        lengthCATG=matrix(),
        scoreCATG=matrix()),
    )


setClass("Roche454QualityFilterMetrics",
    representation=list(
        dataTotal="data.frame",
        dataTCAG="data.frame",
        dataCATG="data.frame"),
    prototype=prototype(
        dataTotal=data.frame(),
        dataTCAG=data.frame(),
        dataCATG=data.frame())
    )


setClass("AnnotatedVariants",
    representation=list(
        annotatedVariants="list"),
    prototype=prototype(
        annotatedVariants=list())
    )


setClass("Breakpoints",
    representation=list(
        seqsC1="list",
	seqsC2="list",
        commonBpsC1="list",
	commonBpsC2="list",
	commonAlignC1="list",
	commonAlignC2="list",
	alignedReadsC1="list",
	alignedReadsC2="list"),
    prototype=prototype(
        seqsC1=list(),
        seqsC2=list(),
        commonBpsC1=list(),
        commonBpsC2=list(),
        commonAlignC1=list(),
        commonAlignC2=list(),
        alignedReadsC1=list(),
        alignedReadsC2=list())
    )



# AVASet: class definition for the AVASet subclass of the original
# Biobase eSet
# extend object by two slots for amplicon data, two slots for variant filtering
# and a slot for the reference sequences
setClass("AVASet",
    contains = "eSet",
    representation=list(
        assayDataAmp="AssayData",
        featureDataAmp="AnnotatedDataFrame",
        referenceSequences="AlignedRead",
        variantFilterPerc="numeric",
	variantFilter="character",
	dirs="character"), 	
    prototype=prototype(
        assayDataAmp=assayDataNew(), 
	featureDataAmp=new("AnnotatedDataFrame"),
	referenceSequences=new("AlignedRead"),
	variantFilterPerc=0,
	variantFilter="",
	dirs="",
	new("VersionedBiobase", versions=c(classVersion("eSet"),
            AVASet="0.1"))
    )
)

# MapperSet: class definition for the MapperSet subclass of the original
# Biobase eSet
# extend object by two slots for variant filtering
setClass("MapperSet",
    contains = "eSet",
    representation=list(
        variantFilterPerc="numeric",
	variantFilter="character", 	
	dirs="character"), 	
    prototype=prototype(
	variantFilterPerc=0,
	variantFilter="",
	dirs="",
	new("VersionedBiobase", versions=c(classVersion("eSet"),
            MapperSet="0.1"))
    )
)

# SFFRead: representation of a single read from a Roche Standard Flowgram Format (SFF) file
setClass("SFFRead",
  representation=list(
    name               = "character",
    flowgramFormat     = "numeric",
    flowChars          = "character",
    keySequence        = "character",
    clipQualityLeft    = "numeric",
    clipQualityRight   = "numeric",
    clipAdapterLeft    = "numeric",
    clipAdapterRight   = "numeric",
    flowgram           = "numeric",
    flowIndexes        = "numeric",
    read               = "DNAString",
    quality            = "BString"
  ),
  prototype=prototype(
    name               = "",
    flowgramFormat     = 1,
    flowChars          = "",
    keySequence        = "",
    clipQualityLeft    = 0,
    clipQualityRight   = 0,
    clipAdapterLeft    = 0,
    clipAdapterRight   = 0,
    flowgram           = 0,
    flowIndexes        = 0,
    read               = DNAString(),
    quality            = BString()
  ),
  validity=function(object) {
    if(length(read(object)) != length(quality(object)) | length(read(object)) != length(flowIndexes(object)))
      return("Read, quality and flowIndexes must be of equal length.")
    if(nchar(flowChars(object)) != length(flowgram(object)))
      stop("Flowgram and flowChars must be of equal length.")
    return(TRUE)
  }
)

# SFFContainer: representation of a complete Roche Standard Flowgram Format (SFF) file
setClass("SFFContainer",
  representation=list(
    name               = "character",
    flowgramFormat     = "numeric",
    flowChars          = "character",
    keySequence        = "character",
    clipQualityLeft    = "numeric",
    clipQualityRight   = "numeric",
    clipAdapterLeft    = "numeric",
    clipAdapterRight   = "numeric", 
    flowgrams          = "list",
    flowIndexes        = "list",
    reads              = "QualityScaledDNAStringSet"
  ),
  prototype=prototype(
    name               = "",
    flowgramFormat     = 1,
    flowChars          = "",
    keySequence        = "",
    clipQualityLeft    = numeric(0),
    clipQualityRight   = numeric(0),
    clipAdapterLeft    = numeric(0),
    clipAdapterRight   = numeric(0),
    flowgrams          = list(),
    flowIndexes        = list(),
    reads              = QualityScaledDNAStringSet(DNAStringSet(), PhredQuality(""))
  ),
  validity=function(object) {
    if(length(reads(object)) != length(flowgrams(object)) |
       length(flowgrams(object)) != length(flowIndexes(sff)))
      return("Number of reads, flowgram length and flow indexes length must be equal")
    return(TRUE)
  }
)

