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
  readLengthHist(object, cutoff=0.99, name)
}

\arguments{
  \item{object}{An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.}
  \item{cutoff}{Reads longer than the cutoff-percent quantile are omitted from the plot.}
  \item{name}{Optional character defining the name to be displayed in the histogram title. Defaults
    to the filename for \link{SFFContainer} and left empty for \link{DNAStringSet} and 
    \link{ShortRead}.}
}

\author{
  Christian Ruckert
}

\keyword{methods}
