\name{alignShortReads}

\alias{alignShortReads}
\alias{alignShortReads,AVASet,BSgenome,character,logical-method}
\alias{alignShortReads,AVASet,BSgenome,character,missing-method}
\alias{alignShortReads,AVASet,BSgenome,missing,logical-method}
\alias{alignShortReads,AVASet,BSgenome,missing,missing-method}
\alias{alignShortReads,DNAStringSet,BSgenome,character,logical-method}
\alias{alignShortReads,DNAStringSet,BSgenome,character,missing-method}
\alias{alignShortReads,DNAStringSet,BSgenome,missing,logical-method}
\alias{alignShortReads,DNAStringSet,BSgenome,missing,missing-method}

\title{Exact alignment of DNA sequences against a reference}

\description{This method aligns given sequences against a given reference
genome using the \code{matchPDict} method. Only exact (no errors) and unique
matches are returned.}

\usage{alignShortReads(object, bsGenome, seqNames, ensemblNotation)}

\arguments{
  \item{object}{The reads that should be aligned agiven either as a
\code{DNAStringSet} or a \code{AVASet} instance. In the latter case
the reference sequences are extracted and aligned.}
  \item{bsGenome}{A \code{bsGenome} instance providing the reference
sequences.}
  \item{seqNames}{The names of the sequences in \code{bsGenome} that should
be used. If omitted, all reference sequences are used.}
  \item{ensemblNotation}{If set to TRUE, \dQuote{chr} is removed from the
reference sequences' names in the returned alignment. Default value is FALSE.}
}

\value{An object of class \code{AlignedRead} or a \code{AVASet} instance.}

\details{All reads are aligned against the reference and its reverse
complement. If the reads are not in 5' to 3' orientation, they should
be reversed before. Note that only exact and unique alignments are
reported. Use \code{matchPDict} directly for more flexibility.}


\seealso{\code{\link{matchPDict}}, \code{\link{DNAStringSet}}, 
\code{\link{AlignedRead}}, \code{\link{AVASet}}
}

\examples{
library("BSgenome.Scerevisiae.UCSC.sacCer2")
reads = DNAStringSet(c(
  "CCGTTCAAAGAGCCCTTGGCCCATAATCCACCGGTT",
  "ATCCTGCCACAGGAGTCCATGGAGGTTTCGCCA"))
alignShortReads(reads, Scerevisiae, seqNames="chrIII")
}

\author{Hans-Ulrich Klein}
\keyword{alignShortReads}


