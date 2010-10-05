\name{sequenceCaptureLinkers}

\alias{sequenceCaptureLinkers}
\alias{sequenceCaptureLinkers,character-method}
\alias{sequenceCaptureLinkers,missing-method}

\title{Retrieve NimbleGen's sequence capture linkers}

\description{This method returns the NimbleGen's linker sequences 
             used with their sequence capture arrays. See pp.29-30
             in the NimbleGen Arrays User's Guide.}

\usage{\S4method{sequenceCaptureLinkers}{character}(name)
       \S4method{sequenceCaptureLinkers}{missing}()}

\arguments{\item{name}{Character vector with linker sequences' names}}

\value{\code{sequenceCaptureLinkers} returns a \code{DNAStringSet} with the
       requested linker sequences.}
     
\details{If the argument \code{name} is omitted, both linker
         sequences are returned.}
     
\seealso{\code{\link{removeLinker}}}

\examples{
    sequenceCaptureLinkers()}

\author{Hans-Ulrich Klein}
\keyword{linker}
