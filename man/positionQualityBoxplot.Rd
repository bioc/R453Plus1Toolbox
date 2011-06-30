\name{positionQualityBoxplot}

\docType{methods}

\alias{positionQualityBoxplot}
\alias{positionQualityBoxplot,QualityScaledDNAStringSet-method}
\alias{positionQualityBoxplot,ShortReadQ-method}
\alias{positionQualityBoxplot,SFFContainer-method}

\title{Boxplot Of The Quality For Each Position}

\description{
  Creates a boxplot of the quality scores over all sequences at each position.
}

\usage{
  positionQualityBoxplot(object, range, binsize=10, 
    xlab=paste("Read position in bp (Bin size: ", binsize, "bp)", sep=""), ylab="Quality score",
    col="firebrick1", ...)
}

\arguments{
  \item{object}{An object of class \link{QualityScaledDNAStringSet}, \link{ShortReadQ} or 
    \link{SFFContainer}.}
  \item{range}{A numeric vector of length one or two. If length one only bases from the first until 
    this position are plotted. If two all bases between these two positions are plotted.}
  \item{binsize}{Number of positions to summarize in one box in the plot.}
  \item{xlab}{The X axis label.}
  \item{ylab}{The Y axis label.}
  \item{col}{The plotting color.}
  \item{\dots}{Arguments to be passed to methods, such as graphical parameters (see \sQuote{par}).}
}

\author{
  Christian Ruckert
}

\keyword{methods}
