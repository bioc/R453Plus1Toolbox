\name{Breakpoints-class}
\Rdversion{1.1}
\docType{class}
\alias{Breakpoints-class}
\alias{alignedReadsC1<-}
\alias{alignedReadsC1<-,Breakpoints,list-method}
\alias{alignedReadsC1}
\alias{alignedReadsC1,Breakpoints-method}
\alias{alignedReadsC2<-}
\alias{alignedReadsC2<-,Breakpoints,list-method}
\alias{alignedReadsC2}
\alias{alignedReadsC2,Breakpoints-method}
\alias{[,Breakpoints,ANY,ANY-method}
\alias{[,Breakpoints,ANY,ANY,ANY-method}
\alias{commonAlignC1<-}
\alias{commonAlignC1<-,Breakpoints,list-method}
\alias{commonAlignC1}
\alias{commonAlignC1,Breakpoints-method}
\alias{commonAlignC2<-}
\alias{commonAlignC2<-,Breakpoints,list-method}
\alias{commonAlignC2}
\alias{commonAlignC2,Breakpoints-method}
\alias{commonBpsC1<-}
\alias{commonBpsC1<-,Breakpoints,list-method}
\alias{commonBpsC1}
\alias{commonBpsC1,Breakpoints-method}
\alias{commonBpsC2<-}
\alias{commonBpsC2<-,Breakpoints,list-method}
\alias{commonBpsC2}
\alias{commonBpsC2,Breakpoints-method}
\alias{length,Breakpoints-method}
\alias{mergeBreakpoints,Breakpoints,missing,list-method}
\alias{mergeBreakpoints,Breakpoints,missing,missing-method}
\alias{mergeBreakpoints,Breakpoints,numeric,missing-method}
\alias{names<-,Breakpoints,ANY-method}
\alias{names,Breakpoints-method}
\alias{plotChimericReads,Breakpoints-method}
\alias{seqsC1<-}
\alias{seqsC1<-,Breakpoints,list-method}
\alias{seqsC1}
\alias{seqsC1,Breakpoints-method}
\alias{seqsC2<-}
\alias{seqsC2<-,Breakpoints,list-method}
\alias{seqsC2}
\alias{seqsC2,Breakpoints-method}

\title{Class "Breakpoints"}
\description{
Container to store chimeric reads that were clustered to putative
breakpoints indicating structural variants. Related information like
breakpoint position or alignment information about the chimeric reads
is stored as well.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Breakpoints", ...)}.
Usually, objects will be created by calling the
\code{\link{detectBreakpoints}} method. It is not intended that users
create objects of this class manually.

All slots of this class can be found twice. One slot name ends with
\dQuote{C1} and the other \dQuote{C2}. The slots labeled with
\dQuote{C2} are empty until \code{\link{mergeBreakpoints}} has been
called and contain information about putativly associated breakpoints
detected by \code{\link{mergeBreakpoints}}.
}
\section{Slots}{
  \describe{
    \item{\code{seqsC1}:}{Object of class \code{"list"} with one data
      frame for each breakpoint. The data frame stores all chimeric reads
      covering the first breakpoint together with the alignment
      information.}
    \item{\code{seqsC2}:}{Object of class \code{"list"} with one data
      frame for each breakpoint. The data frame stores all chimeric reads
      covering the second breakpoint together with the alignment
      information.}
    \item{\code{commonBpsC1}:}{Object of class \code{"list"} with one
      data frame for each breakpoint. The data frame stores the
      consensus breakpoint sequence as well as the breakpoint
      coordinates of the first breakpoint.}
    \item{\code{commonBpsC2}:}{Object of class \code{"list"} with one
      data frame for each breakpoint. The data frame stores the
      consensus breakpoint sequence as well as the breakpoint
      coordinates of the second breakpoint.}
    \item{\code{commonAlignC1}:}{Object of class \code{"list"} with one
      object of class \code{\link{PairwiseAlignmentsSingleSubject-class}}
      for each breakpoint storing the alignments of the chimeric
      reads against the consensus breakpoint sequence for the
      first breakpoint.}
    \item{\code{commonAlignC2}:}{Object of class \code{"list"} with one
      object of class \code{\link{PairwiseAlignmentsSingleSubject-class}}
      for each breakpoint storing the alignments of the chimeric
      reads against the consensus breakpoint sequence for the
      second breakpoint.}
    \item{\code{alignedReadsC1}:}{Object of class \code{"list"} with one
      object of class \code{\link{AlignedRead-class}} storing all
      chimeric reads covering the first breakpoint and their alignments.}
    \item{\code{alignedReadsC2}:}{Object of class \code{"list"} with one
      object of class \code{\link{AlignedRead-class}} storing all
      chimeric reads covering the second breakpoint and their alignments.}
  }
}
\section{Methods}{
  \describe{
    \item{alignedReadsC1<-}{\code{signature(object = "Breakpoints",
	value = "list")}: 
      Setter-method for the alignedReadsC1 slot.}
    \item{alignedReadsC1}{\code{signature(object = "Breakpoints")}:
      Getter-method for the alignedReadsC1 slot.}
    \item{alignedReadsC2<-}{\code{signature(object = "Breakpoints",
	value = "list")}:
      Setter-method for the alignedReadsC2 slot.}
    \item{alignedReadsC2}{\code{signature(object = "Breakpoints")}:
      Getter-method for the alignedReadsC2 slot.}
    \item{commonAlignC1<-}{\code{signature(object = "Breakpoints",
	value =  "list")}: 
      Setter-method for the commonAlignC1 slot.}    
    \item{commonAlignC1}{\code{signature(object = "Breakpoints")}:
      Getter-method for the commonAlignC1 slot.}
    \item{commonAlignC2<-}{\code{signature(object = "Breakpoints",
	value = "list")}:
      Setter-method for the commonAlignC2 slot.}
    \item{commonAlignC2}{\code{signature(object = "Breakpoints")}:
      Getter-method for the commonAlignC2 slot.}
    \item{commonBpsC1<-}{\code{signature(object = "Breakpoints",
	value = "list")}:
      Setter-method for the commonBpsC1 slot.}
    \item{commonBpsC1}{\code{signature(object = "Breakpoints")}:
      Getter-method for the commonBpsC1 slot.}
    \item{commonBpsC2<-}{\code{signature(object = "Breakpoints",
	value = "list")}:
      Setter-method for the commonBpsC2 slot.}
    \item{commonBpsC2}{\code{signature(object = "Breakpoints")}:
      Getter-method for the commonBpsC2 slot.}
    \item{seqsC1<-}{\code{signature(object = "Breakpoints",
	value =	"list")}:
      Setter-method for the seqsC1 slot.}
    \item{seqsC1}{\code{signature(object = "Breakpoints")}:
      Getter-method for the seqsC1 slot.}
    \item{sqsC2<-}{\code{signature(object = "Breakpoints",
	value = "list")}:
      Setter-method for the seqsC2 slot.}
    \item{seqsC2}{\code{signature(object = "Breakpoints")}:
      Getter-method for the seqsC2 slot.}    
    \item{[}{\code{signature(x = "Breakpoints", i = "ANY", j = "ANY")}:
      Subsetting a Breakpoints object.}
    \item{length}{\code{signature(x = "Breakpoints")}:
      Returns the number of breakpoints stored.}
    \item{mergeBreakpoints}{\code{signature(breakpoints = "Breakpoints",
	maxDist = "missing", mergeBPs = "list")}:
      Merge presumably related breakpoints.}
    \item{mergeBreakpoints}{\code{signature(breakpoints = "Breakpoints",
	maxDist = "missing", mergeBPs = "missing")}:
      Merge presumably related breakpoints.}
    \item{mergeBreakpoints}{\code{signature(breakpoints = "Breakpoints",
	maxDist = "numeric", mergeBPs = "missing")}:
      Merge presumably related breakpoints.}
    \item{names<-}{\code{signature(x = "Breakpoints", value = "ANY")}:
      Set the names of the breakpoints.}
    \item{names}{\code{signature(x = "Breakpoints")}:
      Get the names of the breakpoints.}
    \item{plotChimericReads}{\code{signature(brpData = "Breakpoints")}:
      Plot the structural variant and the chimeric reads covering its
      breakpoints.}
    \item{summary}{\code{signature(object = "Breakpoints")}:
      Create a data frame summaring information about all breakpoints.}
    \item{table}{\code{signature(... = "Breakpoints")}:
      Create a frequency table of cluster sizes.}
  }
}

\author{
Hans-Ulrich Klein, Christoph Bartenhagen
}

\seealso{
\code{\link{filterChimericReads}}, \code{\link{detectBreakpoints}},
\code{\link{mergeBreakpoints}}, \code{\link{plotChimericReads}}

}

\keyword{classes}
