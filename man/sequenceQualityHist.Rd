\name{sequenceQualityHist}

\docType{methods}

\alias{sequenceQualityHist}
\alias{sequenceQualityHist,QualityScaledDNAStringSet-method}
\alias{sequenceQualityHist,ShortReadQ-method}
\alias{sequenceQualityHist,SFFContainer-method}

\title{A Histogram Of The Sequence Qualities}

\description{
  This function creates a histogram of the mean qualities of the sequences.
}

\usage{
  sequenceQualityHist(object, name)
}

\arguments{
  \item{object}{An object of class \link{QualityScaledDNAStringSet}, \link{ShortReadQ} or 
    \link{SFFContainer}.}
  \item{name}{Optional character defining the name to be displayed in the histogram title. Defaults
    to the filename for \link{SFFContainer} and left empty for \link{QualityScaledDNAStringSet} and
    \link{ShortReadQ}.}
}

\author{
  Christian Ruckert
}

\keyword{methods}
