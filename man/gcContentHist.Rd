\name{gcContentHist}

\docType{methods}

\alias{gcContentHist}
\alias{gcContentHist,DNAStringSet-method}
\alias{gcContentHist,ShortRead-method}
\alias{gcContentHist,SFFContainer-method}

\title{GC-Content Histogram}

\description{
  This function creates a histogram of the relative GC-content per sequence.
}

\usage{
  gcContentHist(object, name)
}

\arguments{
  \item{object}{An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.}
  \item{name}{Optional character defining the name to be displayed in the histogram title. Defaults
    to the filename for \link{SFFContainer} and left empty for \link{DNAStringSet} and 
    \link{ShortRead}.}
}

\author{
  Christian Ruckert
}

\keyword{methods}
