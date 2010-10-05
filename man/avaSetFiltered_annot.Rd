\name{avaSetFiltered_annot}
\alias{avaSetFiltered_annot}
\docType{data}
\title{
    AVASet variant annotations
}
\description{
    These are example annotations for 4 variants of an AVASet (try \code{data(avaSetFiltered)} to retrieve the corresponding \code{link{AVASet-class}} object). 
    The annotations include affected genes, exons and codons as well as resulting amino acid changes and dbSNP identifiers (if the mutation is already known).
}
\usage{data(avaSetFiltered_annot)}

\format{
    Formal class 'AnnotatedVariants'
}
\source{
    \sQuote{Next-generation sequencing technology reveals a characteristic
    pattern of molecular mutations in 72.8% of chronic myelomonocytic 
    leukemia by detecting frequent alterations in TET2, CBL, RAS, and RUNX1} 
    (Kohlmann A et al., J Clin Oncol. 2010 Aug 20;28(24):3858-65. Epub 2010 Jul 19)
}

\examples{
    data(avaSetFiltered_annot)
}
\keyword{datasets}
