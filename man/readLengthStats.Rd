\name{readLengthStats}

\docType{methods}

\alias{readLengthStats}
\alias{readLengthStats,DNAStringSet-method}
\alias{readLengthStats,ShortRead-method}
\alias{readLengthStats,SFFContainer-method}

\title{Statistics For The Read Lengths}

\description{
  This function returns the mean, median, minimum, maximum and standard deviation of the read 
  lengths over a set of sequences.
}

\usage{
  readLengthStats(object)
}

\arguments{
  \item{object}{
    An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.
  }

}

\value{
  A \link{vector} with five entries: mean, median, min, max and sd.
}

\author{
Christian Ruckert
}

\keyword{methods}
