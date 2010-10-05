\name{getVariantPercentages}

\alias{getVariantPercentages}

\title{Variant coverage}

\description{
    This function computes the coverage for each variant (in forward and/or reverse direction) for all samples. The coverage is defined as the 
    percentual amount of reads that cover a variant.
}

\usage{getVariantPercentages(object, direction="both")}

\arguments{
  \item{object}{An instance of class \code{\link{AVASet-class}} or \code{\link{MapperSet-class}}.}
  \item{direction}{A character indicating the direction ("forward", "reverse" or "both").}
}

\details{
    If the direction was set to "both", the percentages are computed over the sum of both directions. Otherwise it is computed only over the 
    occurences in one direction (forward or reverse). The occurences can be accesses via \code{assayData}.}

\value{
    \code{getVariantPercentages} returns a data frame with all percentages/frequencies for all samples.
}

\seealso{
    \code{\link{setVariantFilter}}.
}

\examples{

# load a (filtered) AVA dataset containing 6 samples, 4 amplicons and 4 variants
data(avaSetFiltered)
avaSetFiltered

# both directions
getVariantPercentages(avaSetFiltered, direction="both")
# this is equivalent to
(assayData(avaSetFiltered)[[1]] + assayData(avaSetFiltered)[[3]]) / (assayData(avaSetFiltered)[[2]] + assayData(avaSetFiltered)[[4]])

# forward direction only
getVariantPercentages(avaSetFiltered, direction="forward")
# this is equivalent to
assayData(avaSetFiltered)[[1]] / assayData(avaSetFiltered)[[2]]

# reverse direction only
getVariantPercentages(avaSetFiltered, direction="reverse")
# this is equivalent to
assayData(avaSetFiltered)[[3]] / assayData(avaSetFiltered)[[4]]


}

\author{Christoph Bartenhagen}
