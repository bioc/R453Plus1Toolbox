\name{setVariantFilter}

\alias{setVariantFilter}

\title{Filters output of variant information}

\description{
    This functions sets the filter to display only those variants, whose amplicon coverage (in percent) in forward and reverse direction in at 
    least one sample is higher than a given value. The coverage is defined as the percentual amount of reads that cover a variant.}

\usage{setVariantFilter(object, filter=0)}

\arguments{
  \item{object}{An instance of an \code{link{AVASet-class}} or \code{\link{MapperSet-class}}.}
  \item{filter}{A filter value between 0 and 1. If two values are given in a vector, the variants are filterd according to the forward (first value) 
	and reverse direction (second value) separately. In this case, a variant has to meet both requirements.}
}

\details{
    Setting the filter affects the assayData and the featureData of the variant slot. See also \code{\link{getVariantPercentages}} for further details.
}

\value{
    \code{setVariantFilter} returns the given \code{link{AVASet-class}}/\code{link{MapperSet-class}} instance with an updated filter value.
}

\seealso{
   \code{link{AVASet-class}},
   \code{link{MapperSet-class}},
   \code{\link{getVariantPercentages}}.
}

\examples{

# load an AVA dataset containing 6 samples, 4 amplicons and 259 variants
data(avaSetExample)
avaSetExample

# use only those variants that are covered by at least 10% of all reads in one sample in both directions together (259 -> 4 variants)
avaSetExample = setVariantFilter(avaSetExample, filter=0.1)
avaSetExample

# use only those variants that are covered by at least 0.1% of all reads in one sample in forward direction
# and by at least 0% in reverse direction (259 -> 6 variants)
avaSetExample = setVariantFilter(avaSetExample, filter=c(0.1, 0))
avaSetExample

# reset filter values to zero
avaSetExample = setVariantFilter(avaSetExample, filter=0)
# or simply
avaSetExample = setVariantFilter(avaSetExample)

}

\author{Christoph Bartenhagen}

