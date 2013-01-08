qualityReportSFF <- function(sfffiles, outfile="qcreport.pdf") {
  basename <- gsub(".pdf", "", outfile, fixed=TRUE)
  tmp <<- environment()
  file.copy(system.file("extdata", "qualityReport.Rnw", package="R453Plus1Toolbox"), paste(basename, ".Rnw", sep=""))
  Sweave(paste(basename, ".Rnw", sep=""))
  texi2dvi(file=paste(basename, ".tex", sep=""), pdf=TRUE)
}
