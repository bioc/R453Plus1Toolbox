\name{flowgramBarplot}

\docType{methods}

\alias{flowgramBarplot}
\alias{flowgramBarplot,SFFRead-method}

\title{Create A Barplot Of The Flow Intensities}

\description{
  This function creates a barplot of the flow intensities with one bar
  for each nucleotide in the flow. With the height giving the measured intensity.
}

\usage{
  flowgramBarplot(x, range=c(0, length(flowgram(x))), xlab="Flow sequence", 
    ylab="Flow intensity", col=c(A="black", C="red", G="blue", T="green"), ...)
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
