\name{MapperSet-class}
\Rdversion{1.1}
\docType{class}
\alias{MapperSet-class}
\alias{annotateVariants,MapperSet-method}
\alias{htmlReport,MapperSet-method}
\alias{setVariantFilter,MapperSet-method}
\alias{getVariantPercentages,MapperSet-method}
\alias{getReadStatus}
\alias{getReadStatus,MapperSet-method}

\title{Class to Contain GS Reference Mapper Output}
\description{	 
    Container to store data imported from a project of Roche's GS Reference Mapper Software.
    It stores all information into a Biobase ExpressionSet.
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{MapperSet(filename)}. 
    While filename is a vector containing all sample directories (i.e. directories that contain the files "mapping/454HCDiffs.txt" and
    "mapping/454NewblerMetrics.txt").
}
\section{Slots}{
  \describe{
    \item{\code{assayData}:}{Object of class \code{AssayData}. Contains the number of reads with the respective difference and the total coverage 
	for every variant in forward and reverse direction.}
    \item{\code{featureData}:}{Object of class \code{AnnotatedDataFrame}. Contains information about the type, location and reference of each
        variant. If available, it shows further Ensembl variant-ids for known SNPs.}
    \item{\code{phenoData}:}{Object of class \code{AnnotatedDataFrame}. By default, the phenoData contains the accession number of every sample.}
    \item{\code{variantFilterPerc}:}{Object of class \code{numeric}. Contains a threshold to display only those variants, whose
        coverage (in percent) in forward and reverse direction in at least one sample is higher than this filter value. See
        \code{\link{setVariantFilter}} for details about setting this value.}
    \item{\code{variantFilter}:}{Object of class \code{character}. Contains a vector of variant names whose
        coverage (in percent) in forward and reverse direction in at least one sample is higher than the filter value(s) in
        \code{variantFilterPerc}.}
    \item{\code{dirs}:}{Object of class \code{character}. Based on a directory given at instantiation of the object, it contains a vector of several 
        directories containing all relevant GS Mapper project files.}
    \item{\code{experimentData}:}{Object of class \code{MIAME}. Contains details of the experiment.}
    \item{\code{annotation}:}{Object of class \code{character} Label associated with the annotation package used in the experiment.}
    \item{\code{protocolData}:}{Object of class \code{AnnotatedDataFrame}. Contains additional information about the samples.}
    \item{\code{.__classVersion__}:}{Object of class \code{Versions}. Remembers the R and R453Toolbox version numbers used to created the
        MapperSet instance.}
  }
}
\section{Extends}{
    Class \code{\linkS4class{eSet}}, directly.
    Class \code{\linkS4class{VersionedBiobase}}, by class "eSet", distance 2.
    Class \code{\linkS4class{Versioned}}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{setVariantFilter(object):}{Sets the filter to display only those variants, whose coverage (in percent) in forward and reverse
        direction in at least one sample is higher than the given value.}
    \item{getVariantPercentages(object)}{Computes the coverage for every variant over all reads (forward and/or reverse) and for each
        sample.}
    \item{annotateVariants(object):}{Annotates given genomic variants. See \code{\link{annotateVariants}} for details.}
    \item{htmlReport(object):}{Exports all (filtered) variant data into a html report. See \code{\link{htmlReport}} for details}
    \item{getReadStatus(object):}{Reads the file "454ReadStatus.txt" in the GSM project directory which contains information about the 
        alignment of each read (chr, pos, strand, etc.)and returns its content in a dataframe.}
	 }
}
\author{Christoph Bartenhagen}

\seealso{
    \code{\link{AVASet}},
    \code{\link{annotateVariants}},
    \code{\link{htmlReport}},
    \code{\link{setVariantFilter}},
    \code{\link{getVariantPercentages}}

}

\examples{

# sum up class structure
showClass("MapperSet")

# load a GS Mapper dataset containing 3 samples and 111 variants
data(mapperSetExample)
mapperSetExample

# show contents of assay, feature and pheno data
assayData(mapperSetExample)
fData(mapperSetExample)
pData(mapperSetExample)
}

\keyword{classes}
