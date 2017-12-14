\name{annotateVariants}

\alias{annotateVariants}
\alias{annotateVariants,AVASet,missing-method}
\alias{annotateVariants,MapperSet,BSgenome-method}
\alias{annotateVariants,MapperSet,missing-method}
\alias{annotateVariants,data.frame,missing-method}

\title{Adds genomic information to variants}

\description{This method annotates given genomic variants (mutations).
Annotation includes affected genes, exons and codons. Resulting amino acid
changes are returned as well as dbSNP identifiers, if the mutation is already
known. All information is fetched from the Ensembl GRCh37 server via biomaRt
using the datasets \code{hsapiens_gene_ensembl} and \code{hsapiens_snp}.}

\usage{annotateVariants(object, bsGenome)}

\arguments{
    \item{object}{A data frame storing variants or an instance of AVASet/MapperSet or a data frame (see details).}
    \item{bsGenome}{An object of class \code{BSGenome} giving the genome to be used as reference sequence
    to calculate amino acid changes. This argument is only applicable when \code{object} is of type 
    \code{\link{MapperSet}}. Default is \sQuote{BSgenome.Hsapiens.UCSC.hg19}. Note that the genome
    should fit  to the Ensembl annotation.}
}

\details{If a data frame is given, the following columns must be
  present:
  \tabular{ll}{
    start \tab genomic start position in the current Ensembl genome\cr
    end \tab genomic end position in the current Ensembl genome\cr
    chromosome \tab chromosome in ensembl notation (i.e. "1", "2", ..., "Y")\cr
    strand \tab "+" or "-" relative to the nucleotide bases given below\cr
    seqRef \tab reference sequence\cr
    seqMut \tab sequence of the observed variant\cr
    seqSur \tab reference sequence extended for 3 bases in both directions
  }
  The rownames of the data frame are used as mutations' names (IDs). See
  examples for a properly defined data drame.
}

\value{An object of class \code{AnnotatedVariants}. Affected genes, transcripts
and exon as well as known SNPs are stored in a list-like structure. See the
documentation of class \code{\link{AnnotatedVariants-class}} for details.}

\seealso{\code{\link{AnnotatedVariants-class}}, \code{\link{AVASet-class}}, 
\code{\link{MapperSet-class}}, \code{\link{htmlReport}}}

\examples{
variants = data.frame(
    start=c(106157528, 106154991,106156184),
    end=c(106157528, 106154994,106156185),
    chromosome=c("4", "4", "4"),
    strand=c("+", "+", "+"),
    seqRef=c("A", "ATAG", "---"),
    seqMut=c("G", "----", "ATA"),
    seqSur=c("TACAGAA", "TTTATAGATA", "AGC---TCC"),
    stringsAsFactors=FALSE)
rownames(variants) = c("snp", "del", "ins")
\dontrun{annotateVariants(variants)}
}

\author{Hans-Ulrich Klein}
\keyword{annotateVariants}
