\name{sff2fastq}

\docType{methods}

\alias{sff2fastq}
\alias{sff2fastq,SFFContainer-method}

\title{Write A SFFContainer Object To A FASTQ File}

\description{
  This function takes a \link{SFFContainer} object and writes it to a file in FASTQ format.
}

\usage{
  sff2fastq(x, outdir, fname)
}

\arguments{
  \item{x}{An object of class \link{SFFContainer}.}
  \item{outdir}{The directory where the file should be stored, defaults to the current working 
    directory.}
  \item{fname}{The name of the file to write. Defaults to the filename slot of the 
    \link{SFFContainer}, with .sff substituted with .fastq.}
}

\author{
  Christian Ruckert
}

\keyword{methods}
