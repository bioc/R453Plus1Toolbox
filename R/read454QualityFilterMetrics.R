.read454QualityFilterMetrics_character <- function(object){

    qfm = read.table(file=object, sep=",",
        colClasses=c("character","character"),
        fill=TRUE, col.names=c("V1", "V2"))

    # Extract the user specified analysis name, i.e. the last part of the
    # analysis directory name.
    ind = qfm[,1] == "analysisName"
    analysisName = unlist(strsplit(qfm[ind, 2], "_"))
    analysisName = analysisName[length(analysisName)]

    # Extract indices
    regInd = qfm[,1] == "region"
    totalRWInd = qfm[,1] == "totalRawWells"
    keySeqInd = qfm[,1] == "keySequence"
    numKeyPassInd = qfm[,1] == "numKeyPass"
    numDotFailedInd = qfm[,1] == "numDotFailed"
    numMixedFailedInd = qfm[,1] == "numMixedFailed"
    numTTSQInd = qfm[,1] == "numTrimmedTooShortQuality"
    numTTSPInd = qfm[,1] == "numTrimmedTooShortPrimer"
    totalPFInd = qfm[,1] == "totalPassedFiltering"

    # Create data frame for aggregated information (both keys)
    dataTotal = data.frame(
        AnalysisName=analysisName,
        Region=qfm[regInd, 2],
        RawWells=as.numeric(qfm[totalRWInd, 2]))
    rownames(dataTotal) = paste(dataTotal$AnalysisName,
        dataTotal$Region, sep="_")

    # Create data frame for TCAG and CATG keys
    dataTmp = data.frame(
        AnalysisName=analysisName,
        Region=as.vector(rbind(qfm[regInd, 2], qfm[regInd, 2])),
        NumKeyPass=as.numeric(qfm[numKeyPassInd, 2]),
        NumDotFailed=as.numeric(qfm[numDotFailedInd, 2]),
        NumMixedFailed=as.numeric(qfm[numMixedFailedInd, 2]),
        NumTrimmedTooShortQuality=as.numeric(qfm[numTTSQInd, 2]),
        NumTrimmedTooShortPrimer=as.numeric(qfm[numTTSPInd, 2]),
        TotalPassedFiltering=as.numeric(qfm[totalPFInd, 2]))
    dataTmp = split(dataTmp, factor(qfm[keySeqInd, 2]))
    dataTCAG = dataTmp$TCAG
    dataCATG = dataTmp$CATG
    rownames(dataTCAG) = paste(dataTCAG$AnalysisName, dataTCAG$Region, sep="_")
    rownames(dataCATG) = paste(dataCATG$AnalysisName, dataCATG$Region, sep="_")

    qfm = new("Roche454QualityFilterMetrics")
    qfm@dataTotal = dataTotal
    qfm@dataTCAG = dataTCAG
    qfm@dataCATG = dataCATG

    return(qfm)
}

setMethod("read454QualityFilterMetrics", signature="character",
    .read454QualityFilterMetrics_character)
