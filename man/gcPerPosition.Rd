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
  gcPerPosition(object, range=0.95, type="l", col="blue", xlab="Position", ylab="Frequency", ...)
}

\arguments{
  \item{object}{An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.}
  \item{range}{An integer vector of length one or two. If length one only bases up to the percent 
    quantile of read lengths defined by the given value are shown. A vector of length two gives the 
    start and end base between which the GC-content is plotted.}
  \item{type}{The type of the plot (see \sQuote{plot}).}
  \item{col}{The plotting color.}
  \item{xlab}{The X axis label.}
  \item{ylab}{The Y axis label.}
  \item{\dots}{Arguments to be passed to methods, such as graphical parameters (see \sQuote{par}).}
}

\author{
  Christian Ruckert
}

\keyword{methods}
