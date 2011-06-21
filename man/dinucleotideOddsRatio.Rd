\name{dinucleotideOddsRatio}

\docType{methods}

\alias{dinucleotideOddsRatio}
\alias{dinucleotideOddsRatio,DNAStringSet-method}
\alias{dinucleotideOddsRatio,ShortRead-method}
\alias{dinucleotideOddsRatio,SFFContainer-method}

\title{Dinucleotide Odds Ratio}

\description{
  This function calculates the dinucleotide odds ratio for each of the sixtheen possible 
  dinucleotides.
}

\usage{
  dinucleotideOddsRatio(object, name)
}

\arguments{
  \item{object}{An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.}
  \item{name}{Optional character defining the name to be displayed in the histogram title. Defaults
    to the filename for \link{SFFContainer} and left empty for \link{DNAStringSet} and 
    \link{ShortRead}.}
}

\details{
  The dinucleotide odds ratio assigns a value between 0 and 2 to each of the sixtheen possible 
  dinucleotides (AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT). For values below 1
  the dinucleotide is under-represented compared to the randomly expected frequency of this 
  dinucleotide in a sequence of the given length and with the given frequencies of the four 
  nucleotides (A, C, G, T). For values above 1 this dinucleotide is over-represented.
}

\value{
  A matrix with sixtheen columns, one for each dinucleotide, containing the dinucleotide odds ratio 
  values for each sequence in a seperate row.
}

\references{
  Schmieder R. (2011) Quality control and preprocessing of metagenomic datasets.
  \emph{Bioinformatics}, 2011 Mar 15;27(6):863-4. 
}

\author{
  Christian Ruckert
}

\keyword{methods}
