\name{detectBreakpoints}

\alias{detectBreakpoints}
\alias{detectBreakpoints,list-method}

\title{Clustering and consensus breakpoint detection for chimeric reads}
\description{
  Given a set of chimeric reads, this methods computes all putative
  breakpoints. First, chimeric reads are clustered such that all reads
  spanning the same breakpoint form a cluster. Then, a consensus
  breakpoint sequence and breakpoint position is computed for each cluster.
}

\usage{
detectBreakpoints(chimericReads, bpDist=100, minClusterSize=4, removeSoftClips=TRUE, bsGenome)
}

\arguments{
  \item{chimericReads}{A list storing chimeric reads as returned by
    \code{\link{filterChimericReads}}. The list must have the format as
    defined by the \code{\link{scanBam}} method.}
  \item{bpDist}{The maximum distance in base pairs between the
    breakpoints of two chimeric reads at which the reads are merge
    to a cluster.}
  
  \item{minClusterSize}{Cluster whose size is below minClusterSize are
    be excluded from breakpoint detection.}
  \item{removeSoftClips}{If true, soft-clipped bases at the beginning or the
    end of a sequence are removed (see details below).}
  \item{bsGenome}{A \code{bsGenome} instance providing the reference
    sequences. If missing, the library
    \code{BSgenome.Hsapiens.UCSC.hg19} is used by default.}
}


\details{
  This method is usually invoked after calling
  \code{\link{filterChimericReads}} and before calling
  \code{\link{mergeBreakpoints}}. It first forms clusters of chimeric
  reads (reads with exactly two local alignments) that span the same
  breakpoint and than computes a consensus breakpoint sequence for each
  cluster.
  
  To carry out a hierarchical clustering, a measure for the distance
  between two chimeric reads must be defined. If reads span different
  chromosomes, their distance is set to infinity. The strand
  information of the local alignments may also indicate that two chimeric
  reads do not span the same breakpoint even if they span the same
  chromosomes. For example, the first reads has two local alignments on
  the positive strand whereas the second read has one local alignment on
  the positive strand and the other on the negative strand. In this case,
  the distance is set to inifinty, too. Finally, the distance measure
  distinguishes between the two breakpoints (sometimes
  called the pathogenic and the reciproce breakpoint) that originate from
  the same structual variant. The distance between a read from the
  pathogenic and a read from the reciproce breakpoint is infinity so that
  two different clusters will emerge. These two related breakpoints can be
  merge later using the \code{\link{mergeBreakpoints}} method. We observed
  that the breakpoints of these two cases often differ by a few ten or even
  a few hundred basepairs.
  
  If the chromosome and strand information between two reads \eqn{x} and
  \eqn{y} are coherent, the Euclidian distance is used:
  \deqn{d(x,y) = (bp(x, ChrA) - bp(y, ChrA))^2 + (bp(x, ChrB) - bp(y, ChrB))^2}
  where \eqn{bp} gives the coordinates of the breakpoint for the given
  read and chromosome. Hierarchical clustering is applied with complete
  linkage and the dendrogram is cutted at a height of \code{bpDist} to
  obtain the final clusters. The \code{bpDist} argument does usually not
  influence the result, because we observed that reads spanning the same
  breakpoint have very little variation (only a few base pairs) in their
  local alignments due to sequencing errors or due to ambiguity caused by
  same/similar sequence of both chromosome near the breakpoint.
  
  Although the given set of reads may belong to the same chimeric DNA,
  their individual breakpoints may differ in a few base pairs. 
  Furthermore, a single read may have more than one possible breakpoint
  if a (small) part of the read was aligned to both parts. \cr
  The following step determines a consensus breakpoint for each cluster.
  It uses the supplied \code{bsGenome} to construct a
  chimeric reference sequence for all possible breakpoints over all reads within
  each cluster. After the reads were realigned to the chimeric reference
  sequences, the one that yields the highest alignment score is taken to 
  represent best the chimeric DNA and its breakpoints. \cr

  As a preprocessing step, \code{detectBreakpoints} offers to remove soft
  clips occuring after the alignment: \cr
  Some reads may contain soft-clipped bases (e.g. linker sequences) at the beginning of the first part
  of the read or at the end of the second part. By default,
  \code{detectBreakpoints} removes these unaligned subsequences
  and adjusts the cigar string, the sequence, the sequence width
  (qwidth) and the local start/end coordinates.
  
}

\value{
  \code{detectBreakpoints} returns an object of class
  \code{breakpoints}, which is a list of breakpoint clusters, which gives
  access to all alignments and consensus breakpoints:  \cr
  \item{seqs}{This \code{\link{IRanges}} \code{\link{DataFrame}} is
    mainly a rearranged version of the alignment input in
    \code{chimericReads}. In addition, it shows the corresponding
    breakpoints and local alignment coordinates.}
  \item{commonBps}{A dataframe listing the breakpoints for both parts
    of the chimeric reference, the associated chromosome, strand and the reference
    sequence itself, including positions "localStart"/"localEnd"
    indicating which part of the reference belongs to which breakpoint.}
  \item{commonAlign}{An object of class
    \code{\link{PairwiseAlignmentsSingleSubject}} of the
    Biostrings package that contains the alignment to the
    (best) consensus reference sequence.}
  \item{alignedReads}{On the basis of \code{commonAlign} and \code{commonBps},
    \code{alignedReads} is an instance of class \code{\link{AlignedRead}}
    containing all aligned reads including their associated chromosomes, strands, and
    positions. Since the reference is a chimeric sequence each read has
    two chromosome and two strand entries.}
}

\author{Hans-Ulrich Klein, Christoph Bartenhagen}

\seealso{
  \code{\link{filterChimericReads}}
  \code{\link{mergeBreakpoints}}
  \code{\link{plotChimericReads}}
}

\examples{
# Construct a small example with three chimeric reads
# (=6 local alignments) in bam format as given by
# aligners such as BWA-SW.
# The first two reads originate from the same case but
# from different strands. The third read originate from
# the reciprocal breakpoint.
library("BSgenome.Scerevisiae.UCSC.sacCer2")
bamReads = list()
bamReads[[1]] = list(
    qname=c("seq1", "seq1", "seq2", "seq2", "seq3", "seq3"),
    flag = as.integer(c(0, 0, 16, 16, 0, 0)),
    rname = factor(c("II", "III", "III", "II", "III", "II")),
    strand = factor(c("+", "+", "-", "-", "+", "+")),
    pos = as.integer(c(99951, 200000, 200000, 99951, 199950, 100001)),
    qwidth = as.integer(c(100, 100, 100, 100, 100, 100)),
    cigar = c("50M50S","50S50M","50S50M","50M50S","50M50S", "50S50M"),
    seq = DNAStringSet(c(
        paste(substr(Scerevisiae$chrII, start=99951, stop=100000),
              substr(Scerevisiae$chrIII, start=200000, stop=200049),
              sep=""),
        paste(substr(Scerevisiae$chrII, start=99951, stop=100000),
              substr(Scerevisiae$chrIII, start=200000, stop=200049),
              sep=""),
        paste(substr(Scerevisiae$chrIII, start=200000, stop=200049),
              substr(Scerevisiae$chrII, start=99951, stop=100000),
              sep=""),
        paste(substr(Scerevisiae$chrIII, start=200000, stop=200049),
              substr(Scerevisiae$chrII, start=99951, stop=100000),
              sep=""),
        paste(substr(Scerevisiae$chrIII, start=199950, stop=199999),
              substr(Scerevisiae$chrII, start=100001, stop=100050),
              sep=""),
        paste(substr(Scerevisiae$chrIII, start=199950, stop=199999),
              substr(Scerevisiae$chrII, start=100001, stop=100050),
              sep="")))
)

bps = detectBreakpoints(bamReads, minClusterSize=1, bsGenome=Scerevisiae)
summary(bps)
table(bps)

mergeBreakpoints(bps)
}
