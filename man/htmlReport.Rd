\name{htmlReport}

\alias{htmlReport}

\title{HTML-Report Builder for the AVASet and MapperSet}

\description{
    This function creates a HTML variant and quality report for a given AVASet or MapperSet instance.
}

\usage{htmlReport(object, annot, blocks=c(), transcripts=c(), sampleCols, minMut=3, dir="HTMLReport", title="Summary")}

\arguments{
  \item{object}{An \code{\link{AVASet-class}} or \code{\link{MapperSet-class}} instance.}
  \item{annot}{An instance of class \code{AnnotatedVariants} you get
    by calling \code{\link{annotateVariants}}. 
    If no such argument is supplied the data will be read from the
    Ensembl database automatically for all variants.}
  \item{blocks}{Character vector of block names for each variant. The
    variants will then be structured into several blocks. If no such list is 
    supplied, the report consists of just one (big) table for all variants.}
  \item{transcripts}{Character vector containing Ensembl transcript-IDs
    that order the according entries on the transcript
    pages. Transcripts given in this argument will appear on top of the
    transcript page.}
  \item{sampleCols}{Character vector of column names of the sample data
    (phenoData) of the AVASet/MapperSet object to filter the sample output 
    on the transcript pages. All columns will be listed if no such argument is given.}
  \item{minMut}{If the value of minMut is greater than zero, the report
    lists only variants, whose coverage for at least one sample is
    higher than minMut (percentage between 0 and 100).}
  \item{dir}{Character with the desired output directory. By defaultm
    the directory "HTMLReport" will be created in the current directory.}
  \item{title}{Heading for the first page with the variant information.}
}

\details{
  The report is structured into two (MapperSet) or three (AVASet)
  parts containing variant and quality information:
  \itemize{
    \item The main page sums up given variant information like the name,
      type, reference gene, position (see \code{\link{fData}},
      \code{\link{annotateVariants}}). \cr
      Using the argument \code{blocks}, the main page can be individually 
      structured by assigning a block name to each variant. \cr
      The main page can be further structured by samples. For a
      given AVASet object, every sample links to another short
      quality report showing only the amplicon coverage for this sample.
    \item Every variant on the main page links to a
      page with further details about the affected genes and transcripts
      (e.g. Ensembl gene-IDs, transcript-IDs, codon sequences, changes of
      amino acids (if coding)). 
    \item Only in case of AVASet object: A quality report shows
      the coverage of every amplicon in forward and/or reverse
      direction. Further plots display the coverage by MID and PTP
      (if this information is given in the pheno data of the object).
  }
}

\seealso{
    \code{\link{annotateVariants}}.
}

\examples{
# note: all examples save the report to the directory "htmlReportExample" in your current R working directory

# load a filtered AVA dataset containing 6 samples, 4 amplicons and 4 variants
# and its variant annotations
data("avaSetFiltered")
data("avaSetFiltered_annot")

# create a full report showing all (unfiltered) information
\dontrun{htmlReport(avaSetFiltered, avaSetFiltered_annot, dir="htmlReportExample", title="htmlReport Example", minMut=0)}

# create a report that emphasizes on samples with variants covered by at least 50% of the reads 
\dontrun{htmlReport(avaSetFiltered, avaSetFiltered_annot, dir="htmlReportExample", title="htmlReport Example", minMut=50)}

# create a report that is structured by the reference genes
library("ShortRead")
refs = sapply(fData(avaSetFiltered)$referenceSeq, function(x) 
    subset(pData(alignData(referenceSequences(avaSetFiltered))), pData(alignData(referenceSequences(avaSetFiltered)))$refSeqID == x)$name)
\dontrun{htmlReport(avaSetFiltered, avaSetFiltered_annot, dir="htmlReportExample", title="htmlReport Example", minMut=0, blocks=refs)}

# create a report whose sample information only lists the sample ids
pData(avaSetFiltered)
sampleCols = "SampleID"
\dontrun{htmlReport(avaSetFiltered, avaSetFiltered_annot, dir="htmlReportExample", title="htmlReport Example", minMut=0, sampleCols=sampleCols)}
}

\author{Christoph Bartenhagen, Hans-Ulrich Klein, Christian Ruckert}
