\name{nucleotideCharts}

\docType{methods}

\alias{nucleotideCharts}
\alias{nucleotideCharts,DNAStringSet-method}
\alias{nucleotideCharts,ShortRead-method}
\alias{nucleotideCharts,SFFContainer-method}

\title{Nucleotide Charts}

\description{
  This function plots the relative frequency of the four bases for each position in the sequences.
}

\usage{
  nucleotideCharts(object, range=0.95, linetypes=c(A="l", C="l", G="l", T="l", N="l"), 
    linecols=c(A="black", C="red", G="blue", T="green", N="grey70"), xlab="Position", 
    ylab="Frequency", ...)
}

\arguments{
  \item{object}{An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.}
  \item{range}{An integer vector of length one or two. If length one only bases up to the percent 
    quantile of read lengths defined by the given value are shown. A vector of length two gives the 
    start and end base between which the nucleotide frequencies are plotted.}
  \item{linetypes}{The line types used for the four nucleotids + N, see (see \sQuote{par})}
  \item{linecols}{The colors in which the four nucleotids + N should be plotted.}
  \item{xlab}{The X axis label.}
  \item{ylab}{The Y axis label.}
  \item{\dots}{Arguments to be passed to methods, such as graphical parameters (see \sQuote{par}).}
}

\author{
  Christian Ruckert
}

\keyword{methods}
