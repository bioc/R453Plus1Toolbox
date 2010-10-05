\name{fDataAmp}

\alias{fDataAmp}

\title{Access the amplicon data of an AVASet.}

\description{
    Similar to \code{fData} of the Biobase ExpressionSet, this function returns the feature data of the amplicon slot of an instance of the 
    AVASet.
}

\usage{
    fDataAmp(object) 
}

\arguments{
  \item{object}{An \code{link{AVASet-class}} object.}
}

\value{
    The feature data of the amplicon slot contains the names, primers, start/end positions and reference sequences of all 
    amplicons (see\code{\link{AVASet-class}} for details). It returns a data frame.
}

\seealso{
    \code{\link{featureDataAmp}}, \code{\link{assayDataAmp}},\code{\link{AVASet-class}}
}

\examples{
    # load an AVA dataset containing 6 samples, 4 amplicons and 259 variants
    data(avaSetExample)
    avaSetExample

    # show contents amplicon feature data
    fDataAmp(avaSetExample)
}

\author{Christoph Bartenhagen}
