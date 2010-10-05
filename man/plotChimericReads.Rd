\name{plotChimericReads}

\alias{plotChimericReads}

\title{Plots chimeric reads}

\description{This function plots a given set of aligned chimeric reads along a reference sequence. It plots the breakpoints of translations or 
 	inversions and marks deletions, insertions and mismatches. Optionally, it displays all base pairs in a given region around the breakpoint.}

\usage{plotChimericReads(brpData, geneSymbols=FALSE, plotMut=TRUE, plotBasePairs=FALSE, maxBasePairs=50, legend=FALSE, title="", 
    col=c("red", "green", "black", "orange"))}

\arguments{
  \item{brpData}{A \code{Breakpoints} object containing the consensus breakpoint of all reads and the consensus reference sequence as returned by the methods 
	\code{\link{detectBreakpoints}} and \code{\link{mergeBreakpoints}}. Since only one plot is made, the
	function will only work for objects of class \code{Breakpoints} having length one.}
  \item{geneSymbols}{Boolean value whether to automatically load and plot the gene symbols from the Ensembl database. Additionally,
	\code{geneSymbols} can be a vector of two strings for an own annotation.}
  \item{plotMut}{Boolean value whether to mark deletions, insertions and mismatches.}
  \item{plotBasePairs}{Optionally, \code{plotChimericReads} displays all base pairs in a given region around the breakpoint (see maxBasePairs).}  
  \item{maxBasePairs}{The maximum number of base pairs to be plotted. Only used in conjunction with \code{plotBasePairs=TRUE}.}
  \item{legend}{A logical value (TRUE/FALSE) whether to plot a legend that explains the colouration of the insertions, deletions, mismatches and 
	breakpoints.}
  \item{title}{A title for the plot.}
  \item{col}{A vector of four colours to draw insertions, deletions, mismatches and breakpoints. In this order, the default colours are "red", 
	"green", "black" and "orange" (use \code{colours()} to see a list of possible values).}
}

\details{
  This method is intended to be run after the pipeline for structural
  variant detection. Therefore, see the methods
  \code{\link{filterChimericReads}}, \code{\link{detectBreakpoints}} and
  \code{\link{mergeBreakpoints}} to correctly preprocess your alignment
  before running \code{plotChimericReads}.
}

\note{
    It is recommended to first create and resize the output device (e.g. the plotting window or a pdf file) before plotting. For example, on Unix 
    systems you may try \code{X11(width=w, height=h)} or \code{pdf(file="plotChimericReads.pdf", width=w, height=h)} for some window width w (e.g. 
    w=12) and window height h (e.g. h=6).
}

\seealso{
  \code{\link{Breakpoints-class}},  
  \code{\link{detectBreakpoints}},
  \code{\link{mergeBreakpoints}}
}

\examples{
# load breakpoint data containing twelve chimeric reads describing an inversion in chromosome 16
data("breakpoints")
breakpoints

# standard plot 
# (only arrangement of reads plotted; breakpoints in orange, deletions
# in red, insertions in green and mismatches in black by default)
plotChimericReads(breakpoints)

# plot base pairs in the breakpoint region (+/- 32bp) 
\dontrun{plotChimericReads(breakpoints, plotBasePairs=TRUE, maxBasePairs=32)}

# use custom colours and display a legend:
# deletions="brown", insertions="blue", mismatches="yellow", breakpoints="gray"
plotChimericReads(breakpoints, col=c("brown", "blue", "yellow", "gray"), legend=TRUE)

}

\author{Christoph Bartenhagen}
\keyword{plotChimericReads, detectBreakpoints, mergeBreakpoints, Breakpoints}
