\name{MapperSet}

\alias{MapperSet}
\alias{MapperSet,character-method}

\title{Creating a MapperSet}

\description{
    This function imports a project of Roche's GS Reference Mapper Software.
    It stores all information into an instance of the Biobase ExpressionSet.
}

\usage{MapperSet(dirs, samplenames)}

\arguments{
  \item{dirs}{A character vector containing all sample directories (i.e. directories that contain the files "mapping/454HCDiffs.txt" (required), 
	"mapping/454ReadStatus.txt" (optional), "mapping/454NewblerMetrics.txt"(optional)).}
  \item{samplenames}{A character vector containing samplenames. The order and number of samplenames must be consistent with the filenames to ensure 
	the correctness of the MapperSet. If no samplenames are given, the filenames are used for naming.}
}

\details{An instance of the MapperSet is derived from the Biobase eSet and thus structured into \cr

1. assayData \cr
\describe{
    \item{variantForwCount/variantRevCount:}{Contain the number of reads with the respective difference in forward/reverse direction.} 
    \item{totalForwCount/totalRevCount:}{Contain the total coverage for every variant in forward/reverse direction.}  
}
2. featureData\cr
\describe{ 
    \item{chromosome, start/end:}{Give the location of each variant.}
    \item{referenceBases/variantBase:}{Show the bases changed in each variant.}
    \item{regName:}{The name of the region (gene) where the variant is located.}
    \item{knownSNP:}{Lists Ensembl variant-ids for known SNPs (if any).}
}
3. phenoData \cr
\describe{
    \item{By default, the phenoData contains the accession number of every sample.}{}
}
}

\value{
    An instance of the MapperSet.
}

\seealso{
    \code{\link{AVASet-class}}
}

\examples{
# load a GS Mapper dataset containing 3 samples and 111 variants
data(mapperSetExample)
mapperSetExample
}

\author{Christoph Bartenhagen}
