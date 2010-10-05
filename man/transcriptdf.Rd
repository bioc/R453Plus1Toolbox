\name{transcriptdf}
\alias{transcriptdf}
\docType{data}
\title{
    Ensembl transcript information for variant plotting
}
\description{
    This example contains Ensembl transcript information (here for transcript 'ENST00000380013') as created, used and returned by the function \code{link{plotVariants}}.
    It mainly holds locations of coding areas within the transcript that will then be highlighted in the variant plot.
}
\usage{data(transcriptdf)}
\format{
  A data frame with 9 observations on the following 5 variables.
  \describe{
    \item{\code{ensembl_transcript_id}}{Ensembl id for the transcript}
    \item{\code{rank}}{exon rank}
    \item{\code{cds_start}}{start of a codon}
    \item{\code{cds_end}}{end of a codon}
    \item{\code{cds_length}}{length of a codon}
  }
}
\examples{
    data(transcriptdf)
}
\keyword{datasets}
