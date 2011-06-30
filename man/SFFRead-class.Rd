\name{SFFRead-class}
\Rdversion{1.1}
\docType{class}
\alias{SFFRead-class}
\alias{read<-,SFFRead,QualityScaledDNAStringSet-method}
\alias{read<-}
\alias{read,SFFRead-method}
\alias{read}
\alias{clipAdapterLeft<-,SFFRead,numeric-method}
\alias{clipAdapterLeft,SFFRead-method}
\alias{clipAdapterRight<-,SFFRead,numeric-method}
\alias{clipAdapterRight,SFFRead-method}
\alias{clipQualityLeft<-,SFFRead,numeric-method}
\alias{clipQualityLeft,SFFRead-method}
\alias{clipQualityRight<-,SFFRead,numeric-method}
\alias{clipQualityRight,SFFRead-method}
\alias{flowgram<-,SFFRead,numeric-method}
\alias{flowgram<-}
\alias{flowgram,SFFRead-method}
\alias{flowgram}
\alias{flowIndexes<-,SFFRead,numeric-method}
\alias{flowIndexes,SFFRead-method}
\alias{name<-,SFFRead,character-method}
\alias{name,SFFRead-method}

\title{Class \code{"SFFRead"}}
\description{
This class is a container for a single read from files in Roche's Standard Flowgram Format (SFF).
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SFFRead", ...)}.
Usually, objects will be created by calling the \code{\link{getRead}} method on an object of class
\code{\linkS4class{SFFContainer}}.
}

\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"} representing the name of the read.}
    \item{\code{read}:}{Object of class \code{"QualityScaledDNAStringSet"} containing the basecalled 
      nucleotide sequence of the read together with the quality scores for each of the bases in the 
      sequence, where the values use the standard -log10 probability scale.}
    \item{\code{clipQualityLeft}:}{Object of class \code{"numeric"} representing the position of the 
      first base after the clipping point for an attached quality sequence. If only a combined 
      (quality+adapter) clipping position is computed it should be stored in clipQualityLeft. If no 
      clipping value is computed the field is set to 0. The position values use 1-based indexing.}
    \item{\code{clipQualityRight}:}{Object of class \code{"numeric"} representing the position of 
      the last base before the clipping point for an attached quality sequence. If only a combined 
      (quality+adapter) clipping position is computed it should be stored in clipQualityRight. If no 
      clipping value is computed the field is set to 0. The position values use 1-based indexing.}
    \item{\code{clipAdapterLeft}:}{Object of class \code{"numeric"} representing the position of the 
      first base after the clipping point for an attached adapter sequence. If only a combined 
      (quality+adapter) clipping position is computed it should be stored in clipQualityLeft. If no 
      clipping value is computed the field is set to 0. The position values use 1-based indexing.}
    \item{\code{clipAdapterRight}:}{Object of class \code{"numeric"} representing the position of 
      the last base before the clipping point for an attached adapter sequence. If only a combined 
      (quality+adapter) clipping position is computed it should be stored in clipQualityRight. If no 
      clipping value is computed the field is set to 0. The position values use 1-based indexing.}
    \item{\code{flowgram}:}{Object of class \code{"numeric"} containing the homopolymer stretch 
      estimates for each flow.}
    \item{\code{flowIndexes}:}{Object of class \code{"numeric"} containing the flow positions for 
      each base in the called sequence, i.e. for each base, the position in the flowgram whose 
      estimate resulted in that base being called.}
  }
}

\section{Methods}{
  \describe{
    \item{read<-}{\code{signature(object = "SFFRead", value = "QualityScaledDNAStringSet")}: 
      Setter-method for the read slot.}
    \item{read}{\code{signature(object = "SFFRead")}: 
      Getter-method for the read slot.}
    \item{clipAdapterLeft<-}{\code{signature(object = "SFFRead", value = "numeric")}: 
      Setter-method for the clipAdapterLeft slot.}
    \item{clipAdapterLeft}{\code{signature(object = "SFFRead")}: 
      Getter-method for the clipAdapterLeft slot.}
    \item{clipAdapterRight<-}{\code{signature(object = "SFFRead", value = "numeric")}: 
      Setter-method for the clipAdapterRight slot.}
    \item{clipAdapterRight}{\code{signature(object = "SFFRead")}: 
      Getter-method for the clipAdapterRight slot.}
    \item{clipQualityLeft<-}{\code{signature(object = "SFFRead", value = "numeric")}: 
      Setter-method for the clipQualityLeft slot.}
    \item{clipQualityLeft}{\code{signature(object = "SFFRead")}: 
      Getter-method for the clipQualityLeft slot.}
    \item{clipQualityRight<-}{\code{signature(object = "SFFRead", value = "numeric")}: 
      Setter-method for the clipQualityRight slot.}
    \item{clipQualityRight}{\code{signature(object = "SFFRead")}: 
      Getter-method for the clipQualityRight slot.}
    \item{flowgram<-}{\code{signature(object = "SFFRead", value = "numeric")}: 
      Setter-method for the flowgram slot.}
    \item{flowgram}{\code{signature(object = "SFFRead")}: 
      Getter-method for the flowgram slot.}
    \item{flowIndexes<-}{\code{signature(object = "SFFRead", value = "numeric")}: 
      Setter-method for the flowIndexes slot.}
    \item{flowIndexes}{\code{signature(object = "SFFRead")}: 
      Getter-method for the flowIndexes slot.}
    \item{name<-}{\code{signature(object = "SFFRead", value = "character")}: 
      Setter-method for the name slot.}
    \item{name}{\code{signature(object = "SFFRead")}: 
      Getter-method for the name slot.}
  }
}

\author{
Christian Ruckert
}

\seealso{
\code{\link{readSFF}}, \code{\linkS4class{SFFContainer}}
}

\examples{
showClass("SFFRead")
}

\keyword{classes}
