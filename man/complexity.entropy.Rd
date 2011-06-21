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
  complexity.entropy(object, name)
}

\arguments{
  \item{object}{An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.}
  \item{name}{Optional character defining the name to be displayed in the histogram title. Defaults
    to the filename for \link{SFFContainer} and left empty for \link{DNAStringSet} and 
    \link{ShortRead}.}
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
