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
  gcContentHist(object, xlab="GC content", ylab="Number of sequences", col="firebrick1", 
    breaks=50, ...)
}

\arguments{
  \item{object}{An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.}
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
