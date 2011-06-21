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
  nucleotideCharts(object, range=0.95, name)
}

\arguments{
  \item{object}{An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.}
  \item{range}{An integer vector of length one or two. If length one only bases up to the percent 
    quantile of read lengths defined by the given value are shown. A vector of length two gives the 
    start and end base between which the nucleotide frequencies are plotted.}
  \item{name}{Optional character defining the name to be displayed in the plot title. Defaults
    to the filename for \link{SFFContainer} and left empty for \link{DNAStringSet} and 
    \link{SFFContainer}.}
}

\author{
  Christian Ruckert
}

\keyword{methods}
