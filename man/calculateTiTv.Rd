\name{calculateTiTv}

\alias{calculateTiTv}
\alias{calculateTiTv,AVASet}
\alias{calculateTiTv,MapperSet}


\title{Calculate transition transversion ratio}

\description{When many point mutations are detected, the ration of transitions
to transversions can be used as quality measure to assess the number of false
positives.}

\usage{
  \S4method{calculateTiTv}{AVASet}(object)
  \S4method{calculateTiTv}{MapperSet}(object)
}

\arguments{
  \item{object}{An instance of AVASet or MapperSet storing the detected variants.}
}

\value{A \code{list} with two elements: A substitution matrix summarizing all
observed substitutions and the transition/transversion ratio.}

\details{For more information about the Ti/Tv ratio see 
http://www.broadinstitute.org/gsa/wiki/index.php/QC_Methods}

\author{Hans-Ulrich Klein}

\examples{
data(avaSetExample)
ava = setVariantFilter(avaSetExample, c(0.03, 0.03))
calculateTiTv(ava)
}

\keyword{TiTv,Transition,Transversion}
