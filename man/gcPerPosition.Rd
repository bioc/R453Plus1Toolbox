\name{gcPerPosition}

\docType{methods}

\alias{gcPerPosition}
\alias{gcPerPosition,DNAStringSet-method}
\alias{gcPerPosition,ShortRead-method}
\alias{gcPerPosition,SFFContainer-method}

\title{GC-Content Per Position}

\description{
  This function plots the GC-content frequency per base position.
}

\usage{
  gcPerPosition(object, range=0.95, name)
}

\arguments{
  \item{object}{An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.}
  \item{range}{An integer vector of length one or two. If length one only bases up to the percent 
    quantile of read lengths defined by the given value are shown. A vector of length two gives the 
    start and end base between which the GC-content is plotted.}
  \item{name}{Optional character defining the name to be displayed in the plot title. Defaults
    to the filename for \link{SFFContainer} and left empty for \link{DNAStringSet} and 
    \link{ShortRead}.}
}

\author{
  Christian Ruckert
}

\keyword{methods}
