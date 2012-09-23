\name{ava2vcf}

\docType{methods}

\alias{ava2vcf}
\alias{ava2vcf,AVASet-method}

\title{Convert an AVASet object into a VCF object}

\description{
  Converts all variants in a given AVASet object into a VCF object and writes it to a file in VCF format if filename is given.
}

\usage{
  ava2vcf(object, filename, annot)
}

\arguments{
  \item{object}{An object of class \link{AVASet}.}
  \item{filename}{The name of the VCF file to write in, if ommitted no file is written.}
  \item{annot}{An object of class \code{AnnotatedVariants}. Optional, if given variants are annotated with informations from dbSNP.}
}

\value{An object of class \code{\link{VCF-class}}}

\examples{
  data("avaSetFiltered")
  vcf <- ava2vcf(avaSetFiltered)
}

\seealso{\code{\link{AnnotatedVariants-class}}, \code{\link{AVASet-class}}, \code{\link{VCF-class}}, \code{\link{writeVcf}}}

\author{
  Christian Ruckert
}

\keyword{methods}
