\name{breakpointTable}

\alias{breakpointTable}
\alias{breakpointTable,Breakpoints-method}

\title{Table of breakpoint support}

\description{Counts the number of chimeric reads supporting each breakpoint.}

\usage{breakpointTable(object, ...)}

\arguments{
  \item{object}{An object of class \code{Breakpoints}.}
  \item{\dots}{Additional arguments passed to \code{\link{table}}, such as
        \code{exclude} or \code{useNA}.}
}

\value{
  An object of class \code{table}. The names correspond to the number of
  supporting chimeric reads per breakpoint, and the values give the number
  of breakpoints with that level of support.
}

\seealso{
  \code{\linkS4class{Breakpoints}}
}

\examples{
# load breakpoint data containing twelve chimeric reads describing an inversion in chromosome 16
data("breakpoints")
breakpoints
breakpointTable(breakpoints)
}


\author{Hans-Ulrich Klein}
\keyword{breakpointTable, Breakpoints}
