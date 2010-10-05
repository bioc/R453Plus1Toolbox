\name{genomeSequencerMIDs}

\alias{genomeSequencerMIDs}
\alias{genomeSequencerMIDs,character-method}
\alias{genomeSequencerMIDs,missing-method}

\title{Retrieve GS multiplex sequences}

\description{This method returns the standard multiplex sequences
             used by the Genome Sequence MID library kits.}

\usage{\S4method{genomeSequencerMIDs}{missing}()
       \S4method{genomeSequencerMIDs}{character}(mid)}

\arguments{\item{mid}{Character vector with multiplex sequences' IDs (MIDs)}}

\value{\code{genomeSequencerMIDs} returns a \code{DNAStringSet} with the
       requested multiplex sequences.}
     
\details{If the argument \code{mid} is omitted, all 14 available
         multiplex sequences are returned.}
     
\seealso{\code{\link{demultiplexReads}}}

\examples{
    genomeSequencerMIDs()
    genomeSequencerMIDs(c("MID1", "MID3"))
}

\author{Hans-Ulrich Klein}
\keyword{demultiplex}


