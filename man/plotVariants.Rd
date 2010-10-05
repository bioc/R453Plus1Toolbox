\name{plotVariants}

\alias{plotVariants}
\alias{plotVariants,AnnotatedVariants,character-method}
\alias{plotVariants,AnnotatedVariants,data.frame-method}

\title{Plots variant positions}

\description{This function illustrates the positions and types of all variants within a given transcript.}

\usage{plotVariants(varData, transcript, legend=TRUE, regions=c(), regLabel="(region label)", col=c("grey","white","black","white"), title="")}

\arguments{
  \item{varData}{An instance of class \code{annotatedVariants} you get by calling \code{\link{annotateVariants}}. 
	If no such argument is supplied the data will be read from the Ensembl database automatically for all variants.}
  \item{transcript}{This argument can either be a string containing an
  Ensembl transcript-id or a data frame containing the associated Ensembl
  transcript information. In the latter case, the data frame requires
  the columns "ensembl_transcript_id", "rank", "cds_start", "cds_end"
  and "cds_length". If only the Ensembl transcript-id is passed, the
  function will return an appropriate data frame.}
  \item{legend}{A logical value (TRUE/FALSE) whether to plot a legend that explains the type of mutation.}
  \item{regions}{A vector with transcript positions (start and end for every region) to be highlighted in the plot. 
	It must have the form \code{regions=c(start1, end1, start2, end2, start3, end3 ...)}.}
  \item{regLabel}{If \code{legend=TRUE} and \code{regions} were specified, it is recommended to give a short label that explains the highlighted 
	region(s).}
  \item{col}{A vector specifying four colours for missense point mutations, silent point mutations, nonsense point mutations and deletions (in this order). By default, these mutations are drawn as grey, white, black arrows and grey diamonds respectively.}
  \item{title}{A title for the plot.}
}

\details{The plot shows only coding regions and distinguishes between
  missense, silent and nonsense point mutations and deletions.}

\seealso{
\code{\link{annotateVariants}}.
}

\value{
  If an Ensembl transcript-id is passed as \code{transcript}
  argument, the function will return a data frame containing the Ensembl
  transcript information required for plotting ("ensembl_transcript_id", "rank", "cds_start", "cds_end"
  and "cds_length"). This data frame may prove useful for future plots with
  the same transcript. \cr
  If such a data frame is passed directly, the function returns nothing.
}
\examples{
# one missense, one nonsense point mutation and one deletion
variants = data.frame(
    row.names=c("missense", "deletion", "nonsense"), 
    start=c(106157528, 106157635, 106193892),
    end=c(106157528, 106157635, 106193892),
    chromosome=c("4", "4", "4"),
    strand=c("+", "+", "+"),
    seqRef=c("A", "G", "C"),
    seqMut=c("G", "-", "T"),
    seqSur=c("TACAGAA", "TAAGCAG", "CGGCGAA"),
    stringsAsFactors=FALSE)

# annotate variants with affected genes, exons and codons (may take a minute to finish)
\dontrun{varAnnot = annotateVariants(variants)}

# plot variants vor transcript "ENST00000380013" and mark region within [700, 1000]
\dontrun{plotVariants(varAnnot, transcript="ENST00000380013", legend=TRUE,regions=c(700,1500), regLabel="interesting region", title="'plotVariants' Example")}
}

\author{Christoph Bartenhagen}
