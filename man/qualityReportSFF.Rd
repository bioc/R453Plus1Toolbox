\name{qualityReportSFF}

\alias{qualityReportSFF}

\title{Function To Create A Quality Report In PDF Format}

\description{
  This function takes a character vector consisting of filenames pointing to files in Roche's SFF format as input
  and creates a quality report in PDF format as output.
}

\usage{qualityReportSFF(sfffiles, outfile)}

\arguments{
  \item{sfffiles}{A character vector of the SFF files to read in.}
  \item{outfile}{The name of the PDF report file created. Defaults to \sQuote{qcreport.pdf} in the current directory.}
}

\details{
  The function uses the qualityReport.Rnw file from the extdata directory of the package and \code{Sweave}
  to create a .tex file which is afterwards converted to .pdf format. In the .Rnw file the following quality
  control functions are used: \code{\link{readLengthStats}}, \code{\link{readLengthHist}}, 
  \code{\link{baseQualityStats}}, \code{\link{baseQualityHist}}, \code{\link{sequenceQualityHist}}, 
  \code{\link{positionQualityBoxplot}}, \code{\link{baseFrequency}}, \code{\link{nucleotideCharts}}, 
  \code{\link{gcContent}}, \code{\link{gcPerPosition}}, \code{\link{gcContentHist}}, 
  \code{\link{complexity.dust}}, \code{\link{complexity.entropy}}, \code{\link{dinucleotideOddsRatio}}.
}

\examples{
  \dontrun{
    file <- system.file("extdata", "SFF", "example.sff", package="R453Plus1Toolbox")
    qualityReportSFF(file, "QualityReport.pdf")
  }
}

\author{Christian Ruckert}
