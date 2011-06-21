\name{baseFrequency}

\docType{methods}

\alias{baseFrequency}
\alias{baseFrequency,DNAStringSet-method}
\alias{baseFrequency,ShortRead-method}
\alias{baseFrequency,SFFContainer-method}

\title{Absolute And Relative Frequency Of The Four Bases.}

\description{
  This function returns the absolute and the relative frequency of the four bases (A, C, G, T).
}

\usage{
  baseFrequency(object)
}

\arguments{
  \item{object}{
    An object of class \link{DNAStringSet}, \link{ShortRead} or \link{SFFContainer}.
  }
}

\details{
  This function makes use of the \link{alphabetFrequency} function from package 
  \code{Biostrings}.
}

\value{
  A \link{data.frame} with two columns containing the absolute and relative frequencies 
  respectively and six rows, one for each of the four bases (A, C, G, T), one for other symbols 
  contained in the reads and one summarizing the five aforementioned rows.
}

\author{
  Christian Ruckert
}

\keyword{methods}
