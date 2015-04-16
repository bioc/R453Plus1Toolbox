\name{removeLinker}

\alias{removeLinker}
\alias{removeLinker,XStringSet,DNAString,logical,numeric,numeric-method}
\alias{removeLinker,XStringSet,DNAString,missing,missing,missing-method}

\title{Remove linker sequences located at the start of short reads}

\description{If linkers are attached during sample preparation, it may
             be useful to remove the linkers' sequences after sequencing.
             This method finds and removes linker sequences that are
	     located at the start of the given reads.
}

\usage{
\S4method{removeLinker}{XStringSet,DNAString,logical,numeric,numeric}(reads, linker, removeReadsWithoutLinker, minOverlap, penalty)
}

\arguments{
  \item{reads}{A \code{DNAStringSet} instance that contains reads
  possibly having linkers at their start site}
  \item{linker}{A \code{DNAString} instance with the linker's sequence}
  \item{removeReadsWithoutLinker}{Whether reads without linkers should
  be removed. Default is FALSE}
  \item{minOverlap}{The minimal score that must be achived when aligning
  the linker. Default is length(linker)/2}
  \item{penalty}{The penalty for substitutions or indels. Default is 2}
}

\value{\code{removeLinker} returns a \code{DNAStringSet} with trimmed reads.}

\details{The best alignment of the linker within the start (length of
  linker + 5) of each given sequence is computed. The followong scoring
  schema is used: Each matching bases scores +1. Each substitution or
  indel scores the given penalty argument (default: penalty=2). There
  are no penalties for gaps and the end of the linker (overlap). An
  alignment is considered as match, if the scores is larger of equal to
  minOverlap (default: minOverlap=round(length(linker)/2)). In cases of
  a successful match, the subsequence from position 1 until the end of
  the linker's alignment is removed.  
}
  

\seealso{\code{\link{sequenceCaptureLinkers}}, \code{\link{DNAStringSet}}, 
  \code{\link{pairwiseAlignment}}
}

\examples{
    linker = sequenceCaptureLinkers()[[1]]
    reads = DNAStringSet(c(
        "CTCGAGAATTCTGGATCCTCAAA",
             "GAATTCTGGATCCTCAAA",
        "CTCGAGAAAAAAAAATCCTCAAA"))
    removeLinker(reads, linker)
}

\author{Hans-Ulrich Klein}
\keyword{removeLinker}


