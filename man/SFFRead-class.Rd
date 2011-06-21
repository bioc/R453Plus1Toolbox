\name{SFFRead-class}
\Rdversion{1.1}
\docType{class}
\alias{SFFRead-class}
\alias{bases<-,SFFRead,DNAString-method}
\alias{bases<-}
\alias{bases,SFFRead-method}
\alias{bases}
\alias{clipAdapterLeft<-,SFFRead,numeric-method}
\alias{clipAdapterLeft,SFFRead-method}
\alias{clipAdapterRight<-,SFFRead,numeric-method}
\alias{clipAdapterRight,SFFRead-method}
\alias{clipQualityLeft<-,SFFRead,numeric-method}
\alias{clipQualityLeft,SFFRead-method}
\alias{clipQualityRight<-,SFFRead,numeric-method}
\alias{clipQualityRight,SFFRead-method}
\alias{flowgramValues<-,SFFRead,numeric-method}
\alias{flowgramValues<-}
\alias{flowgramValues,SFFRead-method}
\alias{flowgramValues}
\alias{flowIndexes<-,SFFRead,numeric-method}
\alias{flowIndexes,SFFRead-method}
\alias{readname<-,SFFRead,character-method}
\alias{readname<-}
\alias{readname,SFFRead-method}
\alias{readname}
\alias{qualityScores<-,SFFRead,numeric-method}
\alias{qualityScores,SFFRead-method}

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
    \item{\code{readname}:}{Object of class \code{"character"} representing the name of the read.}
    \item{\code{bases}:}{Object of class \code{"DNAString"} containing the basecalled nucleotide 
      sequence of the read.}
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
    \item{\code{flowgramValues}:}{Object of class \code{"numeric"} containing the homopolymer 
      stretch estimates for each flow.}
    \item{\code{flowIndexes}:}{Object of class \code{"numeric"} containing the flow positions for 
      each base in the called sequence, i.e. for each base, the position in the flowgram whose 
      estimate resulted in that base being called.}
    \item{\code{qualityScores}:}{Object of class \code{"numeric"} containing the quality scores for
      each of the bases in the sequence, where the values use the standard -log10 probability 
      scale.}
  }
}

\section{Methods}{
  \describe{
    \item{bases<-}{\code{signature(object = "SFFRead", value = "DNAString")}: 
      Setter-method for the bases slot.}
    \item{bases}{\code{signature(object = "SFFRead")}: 
      Getter-method for the bases slot.}
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
    \item{flowgramValues<-}{\code{signature(object = "SFFRead", value = "numeric")}: 
      Setter-method for the flowgramValues slot.}
    \item{flowgramValues}{\code{signature(object = "SFFRead")}: 
      Getter-method for the flowgramValues slot.}
    \item{flowIndexes<-}{\code{signature(object = "SFFRead", value = "numeric")}: 
      Setter-method for the flowIndexes slot.}
    \item{flowIndexes}{\code{signature(object = "SFFRead")}: 
      Getter-method for the flowIndexes slot.}
    \item{readname<-}{\code{signature(object = "SFFRead", value = "character")}: 
      Setter-method for the readname slot.}
    \item{readname}{\code{signature(object = "SFFRead")}: 
      Getter-method for the readname slot.}
    \item{qualityScores<-}{\code{signature(object = "SFFRead", value = "numeric")}: 
      Setter-method for the qualityScores slot.}
    \item{qualityScores}{\code{signature(object = "SFFRead")}: 
      Getter-method for the qualityScores slot.}
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
