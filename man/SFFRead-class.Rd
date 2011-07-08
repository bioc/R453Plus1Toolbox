\name{SFFRead-class}
\Rdversion{1.1}
\docType{class}
\alias{SFFRead-class}
\alias{SFFRead}
\alias{read<-,SFFRead,QualityScaledDNAStringSet-method}
\alias{read<-}
\alias{read,SFFRead-method}
\alias{read}
\alias{flowChars<-,SFFRead,character-method}
\alias{flowChars,SFFRead-method}
\alias{flowgramFormat<-,SFFRead,numeric-method}
\alias{flowgramFormat,SFFRead-method}
\alias{keySequence<-,SFFRead,character-method}
\alias{keySequence,SFFRead-method}
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
\alias{quality<-}
\alias{quality<-,SFFRead,BString-method}
\alias{quality,SFFRead-method}

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
    \item{\code{read}:}{Object of class \code{"DNAString"} containing the basecalled 
      nucleotide sequence of the read.}
    \item{\code{flowgramFormat}:}{Object of class \code{"numeric"} representing the format used to 
      encode each of the flowgram values for each read. Currently, only one flowgram format has been 
      adopted and is coded by the value 1.}
    \item{\code{flowChars}:}{Object of class \code{"character"} containing the array of nucleotide 
      bases ('A', 'C', 'G', 'T') that correspond to the nucleotides used for each flow of each 
      read.}
    \item{\code{keySequence}:}{Object of class \code{"character"} representing the nucleotide bases
      of the key sequence used for these reads.}
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
    \item{\code{quality}:}{Object of class \code{"BString"} containing
      the quality scores for each of the bases in the sequence, where
      the values use the standard -log10 probability scale.}
  }
}

\section{Methods}{
  \describe{
    \item{read<-}{\code{signature(object = "SFFRead", value = "DNAString")}: 
      Setter-method for the read slot.}
    \item{read}{\code{signature(object = "SFFRead")}: 
      Getter-method for the read slot.}
    \item{flowChars<-}{\code{signature(object = "SFFContainer", value = "character")}: 
      Setter-method for the flowChars slot.}
    \item{flowChars}{\code{signature(object = "SFFContainer")}: 
      Getter-method for the flowChars slot.}
    \item{flowgramFormat<-}{\code{signature(object = "SFFContainer", value = "numeric")}: 
      Setter-method for the flowgramFormat slot.}
    \item{flowgramFormat}{\code{signature(object = "SFFContainer")}: 
      Getter-method for the flowgramFormat slot.}
    \item{keySequence<-}{\code{signature(object = "SFFContainer", value = "character")}: 
      Setter-method for the keySequence slot.}
    \item{keySequence}{\code{signature(object = "SFFContainer")}: 
      Getter-method for the keySequence slot.}
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
    \item{quality<-}{\code{signature(object = "SFFRead", value = "BString")}: 
      Setter-method for the quality slot.}
    \item{quality}{\code{signature(object = "SFFRead")}: 
      Getter-method for the quality slot.}
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
