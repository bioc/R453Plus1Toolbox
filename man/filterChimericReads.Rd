\name{filterChimericReads}

\alias{filterChimericReads}
\alias{filterChimericReads,list,RangesList,DNAString,missing,missing-method}
\alias{filterChimericReads,list,RangesList,DNAString,numeric,numeric-method}
\alias{filterChimericReads,list,RangesList,missing,missing,missing-method}
\alias{filterChimericReads,list,RangesList,missing,numeric,numeric-method}
\alias{filterChimericReads,list,missing,DNAString,missing,missing-method}
\alias{filterChimericReads,list,missing,DNAString,numeric,numeric-method}
\alias{filterChimericReads,list,missing,missing,missing,missing-method}
\alias{filterChimericReads,list,missing,missing,numeric,numeric-method}


\title{Extract chimeric reads and apply filtering steps to remove artificial
  chimeric reads.}

\description{Chimeric reads may be caused by sequencing a chromosomal
  aberration or by technical issues during sample preparation. This
  method implements several filter steps to remove false chimeric
  reads.}

\usage{
  \S4method{filterChimericReads}{list,RangesList,DNAString,numeric,numeric}(alnReads, targetRegion, linkerSeq, minDist, dupReadDist)
}

\arguments{
  \item{alnReads}{A list storing the aligned reads as produced by the function \link{scanBam}.}
  \item{targetRegion}{A object of class rangesList containing the target region
  of e.g. a used capture array. The parameter may be omitted in case of
  a non targeted sequencing approach.}
  \item{linkerSeq}{A linker sequence that was used during sample
    preparation. It may be omitted.}
  \item{minDist}{The minimum distance between two local alignments (see details), default 1000}
  \item{dupReadDist}{The maximum distance between the 5 prime start
  position of two duplicated reads (see details), default 1.}
}

\value{A \code{list} containing only filtered chimeric reads. The \code{list}
  has the same structure like the given argument
  \code{alnReads}. Additionally, one element \dQuote{log} with logging
  information of each filtering step is added.}

\details{The following filter steps are performed:
  
  1. All chimeric reads with exactly two local alignments are
  extracted. Reads with more than two local alignments are discarded.
  
  2. If the targetRegion argument is given, chimeric reads must have
  one local alignment at least overlapping the the target region. If
  both local alignments are outside the target region, the read is
  discarded.

  3. If the linkerSeq argument is given, all chimeric reads that have
  the linker sequence between their local alignments are removed. When
  searching the linker sequence, 4 mismatches or indels are allowed and
  the linker sequence must not start or end within the first or last ten
  bases of the read. The function searches for the linkerSeq and for it's
  reverse complement.

  4. Two local alignment of a read must have minDist reads between the
  alignments (if both alignment are on the same chromosome). Otherwise,
  the read seems to span a deletion and not a chromosomal aberration and
  is discarded.

  5. Duplicated reads are removed. Two reads are duplicated, if the lie
  on the same strand and have the same 5 prime start position. Due to
  sequencing and alignment errors, the start position may vary for a
  maximum of dupReadDist bases. In case of duplicated reads, only the longest
  read is kept.

  Reads passing all filtering steps are returned in the \code{list}
  structure as given by the \code{alnReads} argument (as derived from
  the \code{scanBam} method). A data frame with information about the
  number of reads that passed each filter is added to the list.
}

\author{Hans-Ulrich Klein}

\seealso{\code{\link{detectBreakpoints}},
  \code{\link{mergeBreakpoints}}, 
  \code{\link{Breakpoints-class}},
  \code{\link{scanBam}}, \code{\link{sequenceCaptureLinkers}}}

\examples{
library(Rsamtools)
bamFile = system.file("extdata", "SVDetection", "bam", "N01.bam", package="R453Plus1Toolbox")
bam = scanBam(bamFile)
data(captureArray)
linker = sequenceCaptureLinkers("gSel3")[[1]]
filterReads = filterChimericReads(bam, targetRegion=captureArray, linkerSeq=linker)
}
    
\keyword{filterChimericReads}
