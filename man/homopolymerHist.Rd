\name{homopolymerHist}

\docType{methods}

\alias{homopolymerHist}
\alias{homopolymerHist,SFFRead-method}

\title{Create A Histogram Of The Homopolymer Stretches}

\description{
  This function creates a histogram for the different lengths of the homopolymer stretches with one bar
  for each nucleotide in the flow. With the height giving the number of this homopolymer stretch in the flowgram.
}

\usage{
  homopolymerHist(x, range=c(0, length(flowgram(x))), xlab="Homopolymer length", 
    ylab="Number of homopolymers", col=c(A="black", C="red", G="blue", T="green"), ...)
}

\arguments{
  \item{x}{An object of class \code{\link{SFFRead}}.}
  \item{range}{Two positions between which the flows should be plotted.}
  \item{xlab}{The X axis label.}
  \item{ylab}{The Y axis label.}
  \item{col}{The colors in which the four nucleotids should be plotted.}
  \item{\dots}{Arguments to be passed to methods, such as graphical parameters (see \sQuote{par}).}
}

\author{
  Christian Ruckert
}

\keyword{methods}