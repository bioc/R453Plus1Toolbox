\name{writeSFF}

\alias{writeSFF}

\title{Function To Write Files In Roche's .sff Format}

\description{
    This function takes an object of class \code{\link{SFFContainer-class}} and writes its contents into
    a file in Roche's Standard Flowgram Format (SFF) with the given filename.
}

\usage{writeSFF(sffContainer, filename)}

\arguments{
  \item{sffContainer}{An \code{\link{SFFContainer-class}} object.}
  \item{filename}{The name of the file to write into.}
}

\seealso{
  \code{\link{readSFF}}, \code{\linkS4class{SFFContainer}}.
}

\examples{
  file <- system.file("extdata", "SFF", "example.sff", package="R453Plus1Toolbox")
  sffContainer <- readSFF(file)
  sffContainer2 <- sffContainer[1:5] 
  \dontrun{writeSFF(sffContainer2, "output.sff")}
}

\author{Christian Ruckert}
