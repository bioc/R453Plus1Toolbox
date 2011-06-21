\name{complexity.dust}

\docType{methods}

\alias{complexity.dust}
\alias{complexity.dust,DNAStringSet-method}
\alias{complexity.dust,ShortRead-method}
\alias{complexity.dust,SFFContainer-method}

\title{Sequence Complexity Using The DUST Algorithm}

\description{
  This function evaluates the sequence complexity using the DUST algorithm.
}

\usage{
  complexity.dust(object, name)
}

\arguments{
  \item{object}{An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.}
  \item{name}{Optional character defining the name to be displayed in the histogram title. Defaults
    to the filename for \link{SFFContainer} and left empty for \link{DNAStringSet} and 
    \link{ShortRead}.}
}

\details{
  The complexity score is based on how often different trinucleotides occur and is scaled between 0 
  and 100. A sequence of homopolymer repeats (e.g. TTTTTTTTTT) has a score of 100, of dinucleotide 
  repeats (e.g. TATATATATA) has a score around 49, and of trinucleotide repeats (e.g. TAGTAGTAG) has 
  a score around 32. Scores above seven can be considered low-complexity.
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
