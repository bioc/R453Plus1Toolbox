\name{getAlignedReads}
\alias{getAlignedReads}
\alias{getAlignedReads,AVASet-method}

\title{
  Import reads from an Amplicon Variant Analyzer project
}
\description{
  For a given AVASet, this function imports all aligned reads belonging to all (or some selected) amplicons of all samples.
}
\usage{
getAlignedReads(object, amplicons, dir)
}

\arguments{
  \item{object}{An instance of the \code{link{AVASet-class}}.}
  \item{amplicons}{An (optional) character vector of amplicon names as
    mentioned in the amplicon feature data (see \code{\link{fDataAmp}}).}
  \item{dir}{Usually, the method tries to retrieve the path to the AVA
    project from the given \code{link{AVASet}} object. However, if it
    fails to find the directory, \code{dir} can be used to set
    the root directory of the AVA project.}
}

\details{
  This function reports all reads for all samples together. If you want to get the reads for some samples individually, try subsetting your 
  AVASet as in the examples below.
}

\value{
  One \code{\link{DNAStringSet}} that contains all aligned reads for all samples (eventually restricted to some given amplicons).
}

\author{Christoph Bartenhagen}

\seealso{
  \code{\link{AVASet-class}},
  \code{\link{fDataAmp}}
}

\examples{

# load an AVA dataset containing 6 samples, 4 amplicons and 259 variants
data(avaSetExample)

# import all reads for amplicon "TET2_E11.04" of the first sample
avaProjectDir = system.file("extdata", "AVASet", package = "R453Plus1Toolbox")
alnReads = getAlignedReads(avaSetExample[, 1], dir=avaProjectDir, amplicons="TET2_E11.04")
show(alnReads)

}
