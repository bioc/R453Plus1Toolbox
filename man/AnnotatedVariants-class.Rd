\name{AnnotatedVariants-class}
\Rdversion{1.1}
\docType{class}
\alias{AnnotatedVariants-class}
\alias{annotatedVariants}
\alias{annotatedVariants,AnnotatedVariants-method}
\alias{annotatedVariants<-,AnnotatedVariants,list-method}
\alias{names,AnnotatedVariants-method}
\alias{names<-,AnnotatedVariants,character-method}

\title{Class "AnnotatedVariants"}
\description{A class for storing annotation about variants. An object of this class
is returned by the method \code{annotateVariants}. The class has not been designed
to be created by users directly. }

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AnnotatedVariants")}. The
method \code{annotateVariants} returns AnnotatedVariants-objects.
}

\section{Slots}{
  \describe{
    \item{\code{annotatedVariants}:}{Object of class \code{"list"} with one entry for
      each variant.}
  }
}

\section{Methods}{
  \describe{
    \item{annotatedVariants}{\code{signature(object = "AnnotatedVariants")}: Get the 
list with variants.}
    \item{annotatedVariants<-}{\code{signature(object = "AnnotatedVariants", value = "list")}:
Set a new list with variants.}
    \item{names}{\code{signature(x = "AnnotatedVariants")}: Get the
names of the with variants.}
    \item{names<-}{\code{signature(x = "AnnotatedVariants", value = "character")}:
Set the names of the variants.}
  }
}

\details{The list encapsulated by this class has one element for each variant. Each
element is a nested list with the elements \code{genes}, \code{transcripts}, \code{exons}
and \code{snps}. All these elements are data frames listing genes, transcripts, exons or
snps respectively that were affected by the variant. Use the example below to explore
the data frames' contents.}

\author{Hans-Ulrich Klein}

\seealso{
\code{\link{annotateVariants}}, \code{\link{htmlReport}}
}

\examples{
variants = data.frame(
    start=c(106157528, 106154991,106156184),
    end=c(106157528, 106154994,106156185),
    chromosome=c("4", "4", "4"),
    strand=c("+", "+", "+"),
    seqRef=c("A", "ATAG", "---"),
    seqMut=c("G", "----", "ATA"),
    seqSur=c("TACAGAA", "TTTATAGATA", "AGC---TCC"),
    stringsAsFactors=FALSE)
rownames(variants) = c("snp", "del", "ins")
\dontrun{av = annotateVariants(variants)}
\dontrun{annotatedVariants(av)[["snp"]]}
}


\keyword{classes}
