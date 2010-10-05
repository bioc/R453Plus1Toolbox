\name{referenceSequences}

\alias{referenceSequences}

\title{Access the reference sequences of an AVASet}

\description{
    This function give access to a slot of an instance of the AVASet storing information about all reference sequences of the 
    amplicons.
}

\usage{referenceSequences(object)}

\arguments{
  \item{object}{An \code{link{AVASet-class}} object.}
}

\value{
    The data is stored in an object of class \code{AlignedRead} and thus gives information about all reference sequences and 
    their position on a chromosome (if \code{\link{alignShortReads}} has been called before).
}

\seealso{
    \code{\link{alignShortReads}}
}

\examples{

    # load an AVA dataset containing 6 samples, 4 amplicons and 259 variants
    data(avaSetExample)

    referenceSequences(avaSetExample)

}

\author{Christoph Bartenhagen}
