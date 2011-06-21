\name{baseQualityHist}

\docType{methods}

\alias{baseQualityHist}
\alias{baseQualityHist,QualityScaledDNAStringSet-method}
\alias{baseQualityHist,ShortReadQ-method}
\alias{baseQualityHist,SFFContainer-method}

\title{Plot A Histogram Of The Base Qualities.}

\description{
  Create a histogram based on the quality of every single base from all sequences.
}

\usage{
  baseQualityHist(object, name)
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
