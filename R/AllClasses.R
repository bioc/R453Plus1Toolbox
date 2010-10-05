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

