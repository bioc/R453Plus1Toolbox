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
  baseQualityHist(object, xlab="Quality score", ylab="Number of bases", col="firebrick1", breaks=40,
     ...)
}

\arguments{
  \item{object}{An object of class \link{QualityScaledDNAStringSet}, \link{ShortReadQ} or 
    \link{SFFContainer}.}
  \item{xlab}{The X axis label.}
  \item{ylab}{The Y axis label.}
  \item{col}{The plotting color.}
  \item{breaks}{The number of breaks in the histogram (see \sQuote{hist}).}
  \item{\dots}{Arguments to be passed to methods, such as graphical parameters (see \sQuote{par}).}
}

\author{
  Christian Ruckert
}

\keyword{methods}
