\name{avaSetFiltered}
\alias{avaSetFiltered}
\docType{data}
\title{
    Amplicon Variant Analyzer data import
}
\description{
    This is an example of an \code{link{AVASet-class}} object containing the output of Roche's Amplicon Variant Analyzer Software.
    It consists of 6 samples, 4 amplicons and 4 variants. The variants were previously filtered according to the amplicon coverage 
    (see \code{\link{setVariantFilter}} for details about filtering an AVASet object).
}
\usage{data(avaSetFiltered)}

\format{
    Formal class 'AVASet'
}
\source{
    \sQuote{Next-generation sequencing technology reveals a characteristic
    pattern of molecular mutations in 72.8% of chronic myelomonocytic 
    leukemia by detecting frequent alterations in TET2, CBL, RAS, and RUNX1} 
    (Kohlmann A et al., J Clin Oncol. 2010 Aug 20;28(24):3858-65. Epub 2010 Jul 19)
}

\examples{
    data(avaSetFiltered)
    avaSetFiltered
}
\keyword{datasets}
