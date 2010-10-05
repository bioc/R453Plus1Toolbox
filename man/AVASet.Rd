\name{AVASet}

\alias{AVASet}
\alias{AVASet,character-method}

\title{Creating an AVASet}

\description{
    This function imports a project of Roche's Amplicon Variant Analyzer Software. 
    It stores all information into an extended version of the Biobase eSet.
}

\usage{AVASet(dirname)}

\arguments{
  \item{dirname}{The path of the project data (i.e. a directory that contains the files and subdirectories
	"Amplicons/ProjectDef/ampliconsProject.txt", "Amplicons/Results/Variants/currentVariantDefs.txt", "Amplicons/Results/Variants", 
	"Amplicons/Results/Align").}
}

\details{
An AVASet object consists of three slots to store data about \cr
1. variants \cr
\describe{
    \item{variantForwCount/variantRevCount:}{Data frames that contain the number of reads with the respective variant in forward/reverse 
	direction.}
    \item{totalForwCount/totalRevCount:}{Data frames that contain the total coverage for every variant location in forward/reverse 
	direction.}
    \item{referenceSeq:}{Gives the identifier of the reference sequence.}
    \item{variantBase/referenceBases:}{The bases changed in each variant.}
    \item{start/end:}{The position of the variant on the reference sequence.}
    \item{canonicalPattern/name:}{Short identifiers of a variant including the position and the bases changed.}
}
2. amplicons \cr
\describe{
    \item{forwCount/revCount:}{Data frames that contain the number of reads for every amplicon and each sample in 
	forward/reverse direction.}
    \item{primer1,primer2:}{The primer sequences for every amplicon.}
    \item{referenceSeqID:}{The identifier of the reference sequence.}
    \item{targetStart/targetEnd:}{The coordinates of the target region.}
}
3. reference sequences \cr
\describe{
    \item{If additional information has been loaded from Ensembl via \code{alignShortReads}, this slot knows about the chromosome, position and 
	the strand of each reference sequence.}{}
} 

The structure of the variant and amplicon data is derived from the Biobase eSet and thus separated into assayData, phenoData and 
featureData. All information about the reference sequences is stored into an object of class \code{AlignedRead}.

The phenoData of the variants lists the sample-IDs and name, annotation and group of the read data for all samples. If available, the pico titer 
plate(s) (PTP) or MID(s) of each sample are shown as well.
}

\value{
    An instance of the AVASet class.
}

\seealso{
    \code{\link{MapperSet-class}}, \code{\link{alignShortReads}}
}

\examples{
# load an AVA dataset containing 6 samples, 4 amplicons and 259 variants
data(avaSetExample)
avaSetExample
}

\author{Christoph Bartenhagen}
