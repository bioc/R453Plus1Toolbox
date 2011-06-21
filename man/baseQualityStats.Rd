\name{baseQualityStats}

\docType{methods}

\alias{baseQualityStats}
\alias{baseQualityStats,QualityScaledDNAStringSet-method}
\alias{baseQualityStats,ShortReadQ-method}
\alias{baseQualityStats,SFFContainer-method}

\title{Statistics Of Base Quality}

\description{
  This function returns mean, minimum, maximum and standard deviation of the base quality scores 
  over all sequences.
}

\usage{
  baseQualityStats(object)
}

\arguments{
  \item{object}{An object of class \link{QualityScaledDNAStringSet}, \link{ShortReadQ} or 
    \link{SFFContainer}.}
}

\value{
  A numeric vector with four slots: mean, min, max, sd.
}

\author{
  Christian Ruckert
}

\keyword{methods}
