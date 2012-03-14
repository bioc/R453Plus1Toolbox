\name{plotVariants}

\alias{plotVariants}
\alias{plotVariants,AnnotatedVariants,character-method}
\alias{plotVariants,data.frame,data.frame-method}

\title{Plots variant positions}

\description{This function illustrates the positions and types of mutations within a given gene and transcript. The plot shows only coding regions (thus, units are amino acids / codons). The coding 
region is further divided into exons labeled with their rank in the transcript.}

\usage{plotVariants(data, gene, transcript, regions, mutationInfo, horiz=FALSE, cex=1, title="", legend=TRUE)}

\arguments{
  \item{data}{This can be either a data.frame or an instance of class
  \code{annotatedVariants} you get by calling
  \code{\link{annotateVariants}}. A data frame requires the columns
  "label", "pos" "mutation" and "color" specifying an annotation for each
  mutation, its position (a single numeric value), the mutation type and (optionally) the color of
  the mutation (e.g. depicting a certain group identity independent from
  the mutation type). The position needs to be given as amino acids / codons.}
  \item{gene}{A string containing the Ensembl id of the gene of interest
  (required).}
  \item{transcript}{A string containing an Ensembl transcript id
  (optional, but recommended). If no Ensembl transcript-id is passed, it
  is chosen automatically and the function will return an appropriate data frame with all annotated
  transcripts for the given gene.}
  \item{regions}{A data frame having columns "name", "start", "end" and
  "color". The plot will highlight these regions with the given colors
  and print their name in the legend.}
  \item{mutationInfo}{A data frame with annotations for the different
  mutations occurring in the data (optional but recommended). It requires
  the columns "mutation", "legend" and "color". The first column must
  list the exact mutation names that occur in the data column
  "mutation". The column "legend" allows for a more detailed name of the
  mutation that will appear in the legend. The color of each mutation
  type is optional and can also be assigned automatically.}
  \item{horiz}{In more comprehensive datasets, more than one
  mutation may be listed for a single position. If  \code{horiz=FALSE}, these
  overlapping mutations will be aligned vertically. If
  \code{horiz=TRUE} they will be aligned horizontally in groups allowing
  a label for every single mutation.}
  \item{cex}{A numeric value \code{> 0} giving the factor by which the
  labels are magnified relative to the default text size. If the width
  of the device is too small for all mutation labels, their size will be
  scaled automatically until it fits.}
  \item{title}{A title for the plot (optional).}
  \item{legend}{A logical value (TRUE/FALSE) whether to plot a legend
  for the mutations types (see \code{mutationInfo}) and the
  highlighted regions (if any).}
}

\details{The plot will show the coding part of the gene and its exons as
x-axis at the bottom. Mutations will be marked and ordered above according to
their genomic position. The axis units are amino acids / codons (hence, all given
genomic positions should be divided by three if necessary). \cr
Passing a data frame to this function allows a much more individual and detailed
annotation of the mutations like labels, colors and user defined
mutation types. \cr
Passing an instance of class \code{\link{annotatedVariants}} is useful for integration into the
R453Plus1Toolbox pipeline and for compatibility to older versions of the
plot. Then, it will only distinguish deletions and missense, nonsense and silent mutations.
}

\note{Depending on the amount of mutations and the size of the gene, the plot
may not fit into the device or the text may become too small. It is
recommended to carefully select the right size of your device
before starting this function to ensure a well scaled and beautiful
plot. \cr
This function requires the package \code{TeachingDemos} to work, which
can be found at CRAN.}

\seealso{
\code{\link{annotateVariants}}.
}

\value{
  The function will return a data frame containing all Ensembl
  transcript information for the given gene ("ensembl_transcript_id", "rank", "cds_start", "cds_end"
  and "cds_length"). This data frame may prove useful for retrieving a
  Ensembl transcript id for future plots.
}
\examples{

# EXAMPLE 1: Working with intances of class annotatedVariants

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

# plot variants for gene TET2 having the Ensembl id "ENSG00000168769"
# when passing no transcript, the largest transcript annotated in the Ensembl database for this gene will be selected automatically
\dontrun{plotVariants(data=varAnnot, gene="ENSG00000168769", title="plotVariants Example", legend=TRUE)}

# EXAMPLE 2: Working with a data frame

# two missense at one position, one nonsense point mutation and one deletion
# it is possible to assign a color to every single mutation independently from its type
variants = data.frame(
    label=c("A>G","A>G(2)","delG","C>T"),
    pos=c(831,831,867,1437),
    mutation=c("M","M","D","N"),
    color=c("black", "black", "green", "red"),
    stringsAsFactors=FALSE
)

# more detailed names for mutation abbreviations can be passed as mutationInfo
# this is useful for the legend, but can also be generated automatically
mutationInfo = data.frame(
    mutation=c("M","D","S","N"),
    legend=c("Missense","Nonsense","Silent","Deletion"),
    stringsAsFactors=FALSE
)

# regions of interest can be highlighted using the regions parameter
regions = data.frame(
    name = c("region1", "region2"),
    start = c(700, 1400),
    end = c(1000, 1900),
    color = c("red", "blue")
)

# using the horiz parameter, multiple mutations occurring at the same place can be either aligned ...
# ... vertically
\dontrun{plotVariants(data=variants, gene="ENSG00000168769", transcript="ENST00000513237", regions=regions, mutationInfo=mutationInfo, horiz=FALSE, title="'plotVariants' Example", legend=TRUE)}

# ... or horizontally in groups
\dontrun{plotVariants(data=variants, gene="ENSG00000168769", transcript="ENST00000513237", regions=regions, mutationInfo=mutationInfo, horiz=TRUE, title="'plotVariants' Example", legend=TRUE)}

}

\author{Christoph Bartenhagen}
