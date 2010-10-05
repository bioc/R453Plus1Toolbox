\name{readsOnTarget}
\alias{readsOnTarget}
\alias{readsOnTarget,list,RangesList-method}

\title{
Check for each read whether it aligns within the given region.
}

\description{
This methods checks (approximately) whether the given reads align within
a given region. 
}

\usage{
  readsOnTarget(alnReads, targetRegion) 
}

\arguments{
  \item{alnReads}{
  A list as returned by \code{scanBam} storing aligned reads.
}
  \item{targetRegion}{
  The target region as a \code{RangesList}. The chromosome names must
  fit to the chromosome names used in the alignment information of the
  given reads.
}
}

\details{
The detailed alignment information given by the CIGAR strings in .bam
files are ignored by the function. Instead, it is assumed that the whole
read alignes to the reference without indels. This is often not true for
longer read (e.g. generated with Roche 454 Sequencing), but saves
computation time. Hence, this method is useful to approximate the number
of reads that align in the target region of a targeted sequencing experiment.
}


\value{A list with one logical vector for each list entry in alnReads. 
The logical vector indicates for each read whether it overlaps with at least
one base from any target region or not.}

\author{Hans-Ulrich Klein}

\seealso{
\code{\link{scanBam}}
}

\examples{
library(Rsamtools)
bamFile = system.file("extdata", "StructuralVariantDetection", "bam", "N01.bam", package="R453Plus1Toolbox")
bam = scanBam(bamFile)
region = RangesList("11"=IRanges(start=118307205, end=118395936))
targetReads = readsOnTarget(bam, region)
sum(targetReads[[1]])
}


\keyword{readsOnTarget}
