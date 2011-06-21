\name{gcContent}

\docType{methods}

\alias{gcContent}
\alias{gcContent,DNAStringSet-method}
\alias{gcContent,ShortRead-method}
\alias{gcContent,SFFContainer-method}

\title{Calculate The Overall GC-Content}

\description{
  This function calculates the GC-content summarized over all sequences.
}

\usage{
  gcContent(object)
}

\arguments{
  \item{object}{
    An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.
  }
}

\details{
  The GC-content is calculated as follows: \cr
    (#G + #C / #G + #C + #A + #T) * 100 \cr
  Where #G is the number of base G in all sequences.
}

\value{
  A numeric vector of length one containing the overall GC-content.
}

\author{
  Christian Ruckert
}

\keyword{methods}

