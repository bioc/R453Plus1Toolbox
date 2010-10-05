\name{assayDataAmp}

\alias{assayDataAmp}

\title{Access the amplicon data of an AVASet.}

\description{
    Similar to \code{assayData} of the Biobase ExpressionSet, this function returns the assay data of the amplicon slot of an instance of 
    the AVASet.
}

\usage{
    assayDataAmp(object) 
}

\arguments{
  \item{object}{An \code{link{AVASet-class}} object.}
}

\value{
    The assay data of the amplicon slot consists of a list of two data frames with the number of forward and reverse reads of all amplicons for each 
    sample (see \code{\link{AVASet-class}} for details).
}

\seealso{
\code{\link{fDataAmp}}, \code{\link{featureDataAmp}}, \code{\link{AVASet-class}}
}

\examples{
    # load an AVA dataset containing 6 samples, 4 amplicons and 259 variants
    data(avaSetExample)

    # show contents of amplicon assay data
    assayDataAmp(avaSetExample)
}
\author{Christoph Bartenhagen}
