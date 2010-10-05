\name{demultiplexReads}

\alias{demultiplexReads}
\alias{demultiplexReads,XStringSet,XStringSet,missing,logical-method}
\alias{demultiplexReads,XStringSet,XStringSet,missing,missing-method}
\alias{demultiplexReads,XStringSet,XStringSet,numeric,logical-method}
\alias{demultiplexReads,XStringSet,XStringSet,numeric,missing-method}


\title{Performs MID/Multiplex filtering}

\description{Roche's Genome Sequencer allows to load two or more samples
             on one region. To allocate sequences to samples, each sample has 
             a unique multiplex sequence. The multiplex sequence should
             be the prefix of all sequences from that sample. This method
             demultiplexes a given set of sequences according to the given
             multiplex sequences (MIDs).}

\usage{
\S4method{demultiplexReads}{XStringSet,XStringSet,numeric,logical}(reads, mids, numMismatches, trim)}

\arguments{
  \item{reads}{A \code{DNAStringSet} instance that contains reads starting 
               with MIDs}
  \item{mids}{A \code{DNAStringSet} instance that contains the MIDs}
  \item{numMismatches}{The maximal number of mismatches allowed, default 2.}
  \item{trim}{Whether the MIDs should be cutted-out, default TRUE}
}

\value{\code{demultiplexReads} returns a list with one \code{DNAStringSet} 
       instance for each MID.}

\details{All given MIDs must have the same length. The algorithm
  computes the number of mismachtes for each MID. The read is assigned
  to the MID with the lowest number of mismatches. If two or more MIDs
  have the same number of mismachtes, or if the number of mismachtes is
  greater than the given argument \code{numMismachtes}, the read is not
  assigned to any MID. The default number of allowed mismatches is 2.}

\seealso{\code{\link{genomeSequencerMIDs}}, \code{\link{DNAStringSet}}}

\examples{
	library(Biostrings)
    mids = genomeSequencerMIDs(c("MID1", "MID2", "MID3"))
    reads = DNAStringSet(c(
        paste(as.character(mids[["MID1"]]), "A", sep=""),
        paste(as.character(mids[["MID1"]]), "AA", sep=""),
        paste(as.character(mids[["MID2"]]), "AAA", sep="")))
    demultiplexReads(reads, mids)
}

\author{Hans-Ulrich Klein}
\keyword{demultiplex}


