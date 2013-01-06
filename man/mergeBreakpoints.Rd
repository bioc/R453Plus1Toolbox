\name{mergeBreakpoints}
\alias{mergeBreakpoints}

\title{
Identify and merge related breakpoints caused by the same variant.
}
\description{
Structural variation like transversions or inversion cause two
breakpoints. In the context of fusion genes, these are called the
pathogenic breakpoint and the reciproce breakpoint. The method
\code{\link{detectBreakpoints}} processes each breakpoint individually
and does explicitly not put reads from the pathogenic and reciproce
breakpoint into the same cluster. Hence, it is usually sensible to call
this methods afterward to search for related pairs of breakpoints to
gain more confidence about the existence of a structural variation.
}

\usage{
mergeBreakpoints(breakpoints, maxDist, mergeBPs)
}

\arguments{
  \item{breakpoints}{
  An object of class \code{Breakpoints} storing the breakpoints that will
  (potentially) be merged.
}
  \item{maxDist}{
  The maximal distance in basepairs at which two breakpoints will be merged. Default value is 1000.
}
  \item{mergeBPs}{
  An optional list of vectors of length two giving the breakpoints that
  should be merged. If this argument is given, the method will not
  search for related breakpoints.
}
}

\details{
If the \code{maxDist} argument is given, the method compares each pair
of breakpoints and checks, whether the two breakpoints may belong to the
same structural variation. In addition to the spanned chromosomes, the
orientation and the strand information of the reads spanning the
breakpoints are also compared for this purpose. If chrosmosome,
orientation and strand information of two breakpoints go well together,
they will be merged, if  the absolute distance of the breakpoints on
chromosome A plus the absolute distance on chromosome B is smaller or
equal to \code{maxDist}. If one breakpoint has more than one potential
mate breakpoint for merging, it will be merged with the first candidate
and a warning message is printed. The default value of \code{maxDist}
is 1000.\\

If the \code{mergeBPs} argument is given, the method will not search for
related breakpoints but simply merge the given breakpoints.
\code{mergeBPs} must be a \code{list} with vectors of
length two that either contain the names of the indices of the
breakpoints that should be merged.\\

The arguments \code{maxDist} and \code{mergeBPs} cannot be given
together. The given \code{Breakpoints} object must not contain
breakpoints that have been merged before.
}

\value{
An object of class \code{Breakpoints} storing merged and unmerged
breakpoints. 
}

\author{
Hans-Ulrich Klein
}

\seealso{
\code{\link{detectBreakpoints}}, \code{\link{Breakpoints-class}},
\code{\link{plotChimericReads}}
}

\examples{
# Load bam file and filter chimeric reads
library(Rsamtools)
bamFile = system.file("extdata", "SVDetection", "bam", "N01.bam", package="R453Plus1Toolbox")
bam = scanBam(bamFile)
data(captureArray)
linker = sequenceCaptureLinkers("gSel3")[[1]]
filterReads = filterChimericReads(bam, targetRegion=captureArray, linkerSeq=linker)

# detect breakpoints of size >= 3
breakpoints = detectBreakpoints(filterReads, minClusterSize=3)
table(breakpoints)
summary(breakpoints)

# merge breakpoints
breakpoints = mergeBreakpoints(breakpoints)
summary(breakpoints)
}

\keyword{"reciproce breakpoint"}
\keyword{"pathogenic breakpoint"}
