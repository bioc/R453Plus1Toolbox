\name{plotAmpliconCoverage}

\alias{plotAmpliconCoverage}
\alias{plotAmpliconCoverage,AVASet,character,logical-method}
\alias{plotAmpliconCoverage,AVASet,character,missing-method}
\alias{plotAmpliconCoverage,AVASet,missing,logical-method}
\alias{plotAmpliconCoverage,AVASet,missing,missing-method}

\title{Creates a plot visualizing the number of reads per amplicon}

\description{A function for visualizing the number of reads per amplicon
or per MID / pico titer plate.}


\usage{
\S4method{plotAmpliconCoverage}{AVASet,character,logical}(avaSet, type="amplicon", bothDirections=TRUE, cex.names=0.8, cex.axis=0.8, las=3, col=c(rgb(217/255, 214/255, 209/255), rgb(173/255, 38/255, 36/255)), ...)}

\arguments{
  \item{avaSet}{An instance of AVASet.}
  \item{type}{A character vector specifying the type of plot.}
  \item{bothDirections}{A logical value determining whether the plot
  sums forward and reverse reads or shows them separately.}
  \item{cex.names}{Font size of the amplicon name labels.}
  \item{cex.axis}{Font size of axes' labels.}
  \item{las}{Orientation of amplicon name labels.}
  \item{col}{Colors used in the plot.}
  \item{\dots}{Arguments to be passed to methods, such as graphical
  parameters (see \sQuote{par}).}
}

\details{
If the argumnet type is \dQuote{amplicon}, the number of reads for each
amplicon are visualized. In case of a AVASet with one sample, a
barplot with one bar for each amplicon is created. In case of more than
one sample, a boxplot with one box for each amplicon is plotted. If type
is \dQuote{mid}, a boxplot with one box for each MID is created. If type 
is \dQuote{ptp}, a boxplot with one box for each pico titer plate is 
created.}

\seealso{\code{\link{AVASet}}}

\examples{
\dontrun{data(avaSetExample)
plotAmpliconCoverage(avaSetExample)
plotAmpliconCoverage(avaSetExample[,1])}
}

\author{Hans-Ulrich Klein}

\keyword{amplicon coverage}
