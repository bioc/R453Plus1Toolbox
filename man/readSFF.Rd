\name{readSFF}

\alias{readSFF}

\title{Function To Read In Roche's .sff Files}

\description{
    This function reads in files in Roche's Standard Flowgram Format (SFF) and store the contents in 
    an \code{\link{SFFContainer-class}} object.
}

\usage{readSFF(files)}

\arguments{
  \item{files}{The name of the .sff file to read in or a character vector of multiple file names or
    the name of a directory containing .sff files.}
}

\value{
  An object or a list of objects of class \code{SFFContainer} storing all the information from the 
  .sff file(s).
}

\seealso{
  \code{\link{writeSFF}}, \code{\linkS4class{SFFContainer}}.
}

\examples{
  file <- system.file("extdata", "SFF", "example.sff", package="R453Plus1Toolbox")
  sffContainer <- readSFF(file)
  sffContainer
}

\author{Christian Ruckert}
