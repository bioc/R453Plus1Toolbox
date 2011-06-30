\name{readLengthHist}

\docType{methods}

\alias{readLengthHist}
\alias{readLengthHist,DNAStringSet-method}
\alias{readLengthHist,ShortRead-method}
\alias{readLengthHist,SFFContainer-method}

\title{Histogram Of The Read Lengths}

\description{
  This function plots a histogram of the read lengths.
}

\usage{
  readLengthHist(object, cutoff=0.99, xlab="Read length", ylab="Number of sequences", 
    col="firebrick1", breaks=100, ...)
}

\arguments{
  \item{object}{An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.}
  \item{cutoff}{Reads longer than the cutoff-percent quantile are omitted from the plot.}
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
