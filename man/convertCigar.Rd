\name{convertCigar}

\alias{convertCigar}
\alias{extendedCIGARToList}
\alias{listToExtendedCIGAR}

\title{Basic functions for CIGAR strings}

\description{
    These are temporary methods, that are likely to be replaced by
    methods from the Rsamtools package in near future.
}

\usage{
    extendedCIGARToList(cigars)
    listToExtendedCIGAR(cigarList)
}

\arguments{
    \item{cigars}{A character vector with CIGAR strings.}
    \item{cigarList}{A list of converted CIGAR strings as produced by \code{\link{extendedCIGARToList}}}
}

