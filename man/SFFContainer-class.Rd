\name{SFFContainer-class}
\Rdversion{1.1}
\docType{class}
\alias{SFFContainer-class}
\alias{SFFContainer}
\alias{addRead,SFFContainer,SFFRead-method}
\alias{addRead}
\alias{getRead,SFFContainer,character-method}
\alias{getRead}
\alias{clipAdapterLeft<-,SFFContainer,numeric-method}
\alias{clipAdapterLeft<-}
\alias{clipAdapterLeft,SFFContainer-method}
\alias{clipAdapterLeft}
\alias{clipAdapterRight<-,SFFContainer,numeric-method}
\alias{clipAdapterRight<-}
\alias{clipAdapterRight,SFFContainer-method}
\alias{clipAdapterRight}
\alias{clipQualityLeft<-,SFFContainer,numeric-method}
\alias{clipQualityLeft<-}
\alias{clipQualityLeft,SFFContainer-method}
\alias{clipQualityLeft}
\alias{clipQualityRight<-,SFFContainer,numeric-method}
\alias{clipQualityRight<-}
\alias{clipQualityRight,SFFContainer-method}
\alias{clipQualityRight}
\alias{filename<-,SFFContainer,character-method}
\alias{filename<-}
\alias{filename,SFFContainer-method}
\alias{filename}
\alias{flowChars<-,SFFContainer,character-method}
\alias{flowChars<-}
\alias{flowChars,SFFContainer-method}
\alias{flowChars}
\alias{flowgramFormat<-,SFFContainer,numeric-method}
\alias{flowgramFormat<-}
\alias{flowgramFormat,SFFContainer-method}
\alias{flowgramFormat}
\alias{flowgrams<-,SFFContainer,list-method}
\alias{flowgrams<-}
\alias{flowgrams,SFFContainer-method}
\alias{flowgrams}
\alias{flowIndexes<-,SFFContainer,list-method}
\alias{flowIndexes<-}
\alias{flowIndexes,SFFContainer-method}
\alias{flowIndexes}
\alias{keySequence<-,SFFContainer,character-method}
\alias{keySequence<-}
\alias{keySequence,SFFContainer-method}
\alias{keySequence}
\alias{reads<-,SFFContainer,QualityScaledDNAStringSet-method}
\alias{reads<-}
\alias{reads,SFFContainer-method}
\alias{reads}
\alias{[,SFFContainer,ANY,ANY-method}

\title{Class \code{"SFFContainer"}}
\description{
This class is a container for data from files in Roche's Standard Flowgram Format (SFF).
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SFFContainer", ...)}.
Usually, objects will be created by calling the \code{\link{readSFF}} method on a file in SFF 
format.
}

\section{Slots}{
  \describe{
    \item{\code{filename}:}{Object of class \code{"character"} containing the name of the file this 
      SFFContainer was created from.}
    \item{\code{flowgramFormat}:}{Object of class \code{"numeric"} representing the format used to 
      encode each of the flowgram values for each read. Currently, only one flowgram format has been 
      adopted and is coded by the value 1.}
    \item{\code{flowChars}:}{Object of class \code{"character"} containing the array of nucleotide 
      bases ('A', 'C', 'G', 'T') that correspond to the nucleotides used for each flow of each 
      read.}
    \item{\code{keySequence}:}{Object of class \code{"character"} representing the nucleotide bases
      of the key sequence used for these reads.}
    \item{\code{clipQualityLeft}:}{Object of class \code{"numeric"} representing the position of the 
      first base after the clipping point for an attached quality sequence for each read. If only a 
      combined (quality+adapter) clipping position is computed it should be stored in 
      clipQualityLeft. If no clipping value is computed the field is set to 0. The position values 
      use 1-based indexing.}
    \item{\code{clipQualityRight}:}{Object of class \code{"numeric"} representing the position of 
      the last base before the clipping point for an attached quality sequence for each read. If 
      only a combined (quality+adapter) clipping position is computed it should be stored in 
      clipQualityRight. If no clipping value is computed the field is set to 0. The position values 
      use 1-based indexing.}
    \item{\code{clipAdapterLeft}:}{Object of class \code{"numeric"} representing the position of the 
      first base after the clipping point for an attached adapter sequence for each read. If only a 
      combined (quality+adapter) clipping position is computed it should be stored in 
      clipQualityLeft. If no clipping value is computed the field is set to 0. The position values 
      use 1-based indexing.}
    \item{\code{clipAdapterRight}:}{Object of class \code{"numeric"} representing the position of 
      the last base before the clipping point for an attached adapter sequence for each read. If 
      only a combined (quality+adapter) clipping position is computed it should be stored in 
      clipQualityRight. If no clipping value is computed the field is set to 0. The position values 
      use 1-based indexing.}
    \item{\code{flowgrams}:}{Object of class \code{"list"} containing the homopolymer stretch 
      estimates for each flow using one list item for each read.}
    \item{\code{flowIndexes}:}{Object of class \code{"list"} containing the flow positions for 
      each base in the called sequence, i.e. for each base, the position in the flowgram whose 
      estimate resulted in that base being called. Each read has its own list item.}
    \item{\code{reads}:}{Object of class \code{"QualityScaledDNAStringSet"} containing the 
      basecalled nucleotide sequences of each read together with the quality scores for each of the 
      bases in the sequence using the standard -log10 probability scale.}
  }
}

\section{Methods}{
  \describe{
    \item{addRead}{\code{signature(object = "SFFContainer", read = "SFFRead")}: Adds an object of 
      class \code{\linkS4class{SFFRead}} to the \code{\linkS4class{SFFContainer}}}
    \item{getRead}{\code{signature(object = "SFFContainer", readname = "character")}: Returns the 
      read with the given name as an object of class \code{\linkS4class{SFFRead}}.}
    \item{clipAdapterLeft<-}{\code{signature(object = "SFFContainer", value = "numeric")}: 
      Setter-method for the clipAdapterLeft slot.}
    \item{clipAdapterLeft}{\code{signature(object = "SFFContainer")}: 
      Getter-method for the clipAdapterLeft slot.}
    \item{clipAdapterRight<-}{\code{signature(object = "SFFContainer", value = "numeric")}: 
      Setter-method for the clipAdapterRight slot.}
    \item{clipAdapterRight}{\code{signature(object = "SFFContainer")}: 
      Getter-method for the clipAdapterRight slot.}
    \item{clipQualityLeft<-}{\code{signature(object = "SFFContainer", value = "numeric")}: 
      Setter-method for the clipQualityLeft slot.}
    \item{clipQualityLeft}{\code{signature(object = "SFFContainer")}: 
      Getter-method for the clipQualityLeft slot.}
    \item{clipQualityRight<-}{\code{signature(object = "SFFContainer", value = "numeric")}: 
      Setter-method for the clipQualityRight slot.}
    \item{clipQualityRight}{\code{signature(object = "SFFContainer")}: 
      Getter-method for the clipQualityRight slot.}
    \item{filename<-}{\code{signature(object = "SFFContainer", value = "character")}: 
      Setter-method for the filename slot.}
    \item{filename}{\code{signature(object = "SFFContainer")}: 
      Getter-method for the filename slot.}
    \item{flowChars<-}{\code{signature(object = "SFFContainer", value = "character")}: 
      Setter-method for the flowChars slot.}
    \item{flowChars}{\code{signature(object = "SFFContainer")}: 
      Getter-method for the flowChars slot.}
    \item{flowgramFormat<-}{\code{signature(object = "SFFContainer", value = "numeric")}: 
      Setter-method for the flowgramFormat slot.}
    \item{flowgramFormat}{\code{signature(object = "SFFContainer")}: 
      Getter-method for the flowgramFormat slot.}
    \item{flowgrams<-}{\code{signature(object = "SFFContainer", value = "list")}: 
      Setter-method for the flowgrams slot.}
    \item{flowgrams}{\code{signature(object = "SFFContainer")}: 
      Getter-method for the flowgrams slot.}
    \item{flowIndexes<-}{\code{signature(object = "SFFContainer", value = "list")}: 
      Setter-method for the flowIndexes slot.}
    \item{flowIndexes}{\code{signature(object = "SFFContainer")}: 
      Getter-method for the flowIndexes slot.}
    \item{keySequence<-}{\code{signature(object = "SFFContainer", value = "character")}: 
      Setter-method for the keySequence slot.}
    \item{keySequence}{\code{signature(object = "SFFContainer")}: 
      Getter-method for the keySequence slot.}
    \item{reads<-}{\code{signature(object = "SFFContainer", value = "QualityScaledDNAStringSet")}: 
      Setter-method for the reads slot.}
    \item{reads}{\code{signature(object = "SFFContainer")}: 
      Getter-method for the reads slot.}
    \item{[}{\code{signature(x = "SFFContainer", i = "ANY", j = "ANY")}:
      Subsetting a SFFContainer object.}
  }
}

\author{
  Christian Ruckert
}

\seealso{
  \code{\link{readSFF}}, \code{\linkS4class{SFFRead}}
}

\examples{
  showClass("SFFContainer")
}

\keyword{classes}
