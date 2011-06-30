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
  sequenceQualityHist(object, xlab="Mean of quality scores per sequence", 
    ylab="Number of sequences", col="firebrick1", ...)
}

\arguments{
  \item{object}{An object of class \link{QualityScaledDNAStringSet}, \link{ShortReadQ} or 
    \link{SFFContainer}.}
  \item{xlab}{The X axis label.}
  \item{ylab}{The Y axis label.}
  \item{col}{The plotting color.}
  \item{\dots}{Arguments to be passed to methods, such as graphical parameters (see \sQuote{par}).}
}

\author{
  Christian Ruckert
}

\keyword{methods}
