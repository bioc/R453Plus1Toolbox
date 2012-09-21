\name{AVASet}

\alias{AVASet}
\alias{AVASet,character,missing,missing,missing,missing,missing,missing-method}
\alias{AVASet,character,character,missing,missing,missing,missing,missing-method}
\alias{AVASet,character,missing,character,character,character,character,character-method}
\alias{AVASet,character,missing,character,character,character,missing,missing-method}


\title{Creating an AVASet}

\description{
    This function imports a project of Roche's Amplicon Variant Analyzer (AVA) Software. 
    It stores all information into an extended version of the Biobase eSet.
}

\usage{AVASet(dirname, avaBin, file_sample, file_amp, file_reference, file_variant, file_variantHits)}

\arguments{
  \item{dirname}{The path of the AVA project. \cr 
	Without AVA-CLI (AVA version < 2.6): A directory that contains the files and subdirectories
	"Amplicons/ProjectDef/ampliconsProject.txt", "Amplicons/Results/Variants/currentVariantDefs.txt", "Amplicons/Results/Variants", 
	"Amplicons/Results/Align". \cr 
	Using AVA-CLI (recommended): Path usually ends with directory "projectfolder"}
  \item{avaBin}{The directory containing the AVA-CLI binary doAmplicon (usually "bin" in the AVA installation directory)}
  \item{file_sample}{Sample information exported with the AVA-CLI. File has to be in CSV format.}
  \item{file_amp}{Amplicons exported with the AVA-CLI. File has to be in CSV format.}
  \item{file_reference}{Reference sequences exported with the AVA-CLI. File has to be in CSV format.}
  \item{file_variant}{Variant information exported with the AVA-CLI. File has to be in CSV format.}
  \item{file_variantHits}{Report of variant hits exported with the AVA-CLI. File has to be in CSV format.}
}

\details{
The five arguments for AVA command line interface (AVA-CLI) exports are optional and useful for exported projects, when no AVA software is installed. 
For exporting, start the AVA-CLI with the command "doAmplicon" and use the commands "open", then "list sample", "list amplicon", "list reference", 
"list variant" and "report variantHits". See \code{\link{AVASet-class}} for more details. \cr
Giving only a project directory and the path to the AVA-CLI binary doAmplicon, \code{AVASet} will 
import all information by accessing the AVA-CLI from within R. 

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
plate(s) (PTP) or MID(s) of each sample are shown as well (using the AVA-CLI, PTPs and MIDs cannot be importet at the moment).
}

\value{
    An instance of the AVASet class.
}
\note{
    It is recommended to use the import via AVA-CLI access. Although deprecated, the import for projects created with older version of the AVA software (< v2.6) is still possible.  
}

\seealso{
    \code{\link{AVASet-class}}, \code{\link{MapperSet-class}}, \code{\link{alignShortReads}}
}

\examples{
# Loading a project from AVA version < 2.6:
# Load an AVA dataset containing 6 samples, 4 amplicons and 259 variants
data(avaSetExample)
avaSetExample

# Loading exported data, that was exported via AVA-CLI
# Load an AVA dataset containing 6 samples, 4 amplicons and 222 variants
# by specifying each file exported from the AVA-CLI
projectDir = system.file("extdata", "AVASet_doAmplicon", package="R453Plus1Toolbox")
avaSetExample = AVASet(dirname=projectDir, file_sample="sample.csv", file_amp="amp.csv", file_reference="reference.csv", file_variant="variant.csv", file_variantHits="variantHits.csv")
avaSetExample

# In case AVA software is installed:
# Saying, for example, the AVA software was installed to the directory "/home/User/AVA",
# the easiest way to import a project via AVA-CLI would look like:
# avaSetExample = AVASet(dirname="myProjectDir", avaBin="/home/User/AVA/bin")


}

\author{Christoph Bartenhagen}
