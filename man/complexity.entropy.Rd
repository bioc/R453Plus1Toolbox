\name{complexity.entropy}

\docType{methods}

\alias{complexity.entropy}
\alias{complexity.entropy,DNAStringSet-method}
\alias{complexity.entropy,ShortRead-method}
\alias{complexity.entropy,SFFContainer-method}

\title{Sequence Complexity Using The Shannon-Wiener Algorithm}

\description{
  This function evaluates the sequence complexity using the Shannon-Wiener Algorithm.
}

\usage{
  complexity.entropy(object, xlab="Complexity score (0=low, 100=high)", ylab="Number of sequences", 
    xlim=c(0, 100), col="firebrick1", breaks=100, ...)
}

\arguments{
  \item{object}{An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.}
  \item{xlab}{The X axis label.}
  \item{ylab}{The Y axis label.}
  \item{xlim}{The limits of the X axis.}
  \item{col}{The plotting color.}
  \item{breaks}{The number of breaks in the histogram (see \sQuote{hist}).}
  \item{\dots}{Arguments to be passed to methods, such as graphical parameters (see \sQuote{par}).}
}

\details{
  The entropy approach evaluates the entropy of trinucleotides in a sequence. The entropy values are 
  scaled from 0 to 100 and lower entropy values imply lower complexity. A sequence of homopolymer 
  repeats (e.g. TTTTTTTTTT) has an entropy value of 0, of dinucleotide repeats (e.g. TATATATATA) has 
  an entropy value around 16, and of trinucleotide repeats (e.g. TAGTAGTAG) has an entropy value 
  around 26. Scores below 70 can be considered low-complexity.
}

\value{
  A numeric vector containing the complexity score for each sequence.
}

\references{
  Schmieder R. (2011) Quality control and preprocessing of metagenomic datasets.
  \emph{Bioinformatics}, 2011 Mar 15;27(6):863-4. 
}

\author{
  Christian Ruckert
}

\keyword{methods}
