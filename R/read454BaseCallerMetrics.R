.read454BaseCallerMetrics_character <- function(object){

    bcm = read.table(file=object, sep=",",
                     colClasses=c("character","character", "character"),
                     fill=TRUE, col.names=c("V1", "V2", "V3"))

    # Extract the user specified analysis name, i.e. the last part of the
    # analysis directory name.
    ind = bcm[,1] == "Rig Result Directory:"
    analysisName = unlist(strsplit(bcm[ind, 2], "_"))
    analysisName = analysisName[length(analysisName)]

    # Extract the pico titer plate ID.
    ind = bcm[,1] == "PTPBarCode:"
    picoTiterPlateID = as.numeric(bcm[ind,2])

    # Extract reads' length and bases' quality for each
    # region and key
    regInd = which(bcm[,1] == "region:")
    keyInd = which(bcm[,1] == "key:")
    lengthInd = which(bcm[,1] == "Length")
    scoreInd = which(bcm[,1] == "Score")

    dataTCAG = data.frame()
    dataCATG = data.frame()
    lengthTCAG = matrix()
    scoreTCAG = matrix()
    lengthCATG = matrix()
    scoreCATG = matrix()
    
    for (r in 1:length(regInd)) {
        l = as.numeric(bcm[(lengthInd[r] + 1):(scoreInd[r]-1), 1])
        cl = as.numeric(bcm[(lengthInd[r] + 1):(scoreInd[r]-1), 2])

        if (r == length(regInd)) {
          endS = nrow(bcm)        
        } else {
          endS = regInd[r+1] - 1
        }
        s = as.numeric(bcm[(scoreInd[r] + 1):endS, 1])
        cs = as.numeric(bcm[(scoreInd[r] + 1):endS, 2])

        if (bcm[keyInd[r], 2] == "TCA") {
            ind = nrow(dataTCAG) + 1

            m = max(l, nrow(lengthTCAG))
            old = lengthTCAG
            lengthTCAG = matrix(0, nrow=m, ncol=ind)
            lengthTCAG[1:nrow(old), ind-1] = old
            lengthTCAG[l, ind] = cl
        
            m = max(s+1, nrow(scoreTCAG))
            old = scoreTCAG
            scoreTCAG = matrix(0, nrow=m, ncol=ind)
            scoreTCAG[1:nrow(old), ind-1] = old
            scoreTCAG[s+1, ind] = cs

            allSeqsLengths = rep(1:nrow(lengthTCAG), lengthTCAG[, ind])
            allBaseScores = rep(1:nrow(scoreTCAG), scoreTCAG[, ind])
            dataTCAG[ind, "AnalysisName"] = analysisName
            dataTCAG[ind, "PicoTiterPlateID"] = picoTiterPlateID
            dataTCAG[ind, "Region"] = bcm[regInd[r], 2]
            dataTCAG[ind, "NumberOfReads"] = length(allSeqsLengths)
            dataTCAG[ind, "NumberOfBases"] = length(allBaseScores)
            dataTCAG[ind, "ReadLengthMean"] = mean(allSeqsLengths)
            dataTCAG[ind, "ReadLengthStd"] = sd(allSeqsLengths)
            dataTCAG[ind, "ReadLengthMedian"] = median(allSeqsLengths)
            dataTCAG[ind, "ReadLengthMAD"] = mad(allSeqsLengths, constant=1)
            dataTCAG[ind, "ReadScoreMean"] = mean(allBaseScores)
            dataTCAG[ind, "ReadScoreStd"] = sd(allBaseScores)
            dataTCAG[ind, "ReadScoreMedian"] = median(allBaseScores)
            dataTCAG[ind, "ReadScoreMAD"] = mad(allBaseScores, constant=1)
            rm(allSeqsLengths)
            rm(allBaseScores)
 
        } else if (bcm[keyInd[r], 2] == "CAT") {
            ind = nrow(dataCATG) + 1

            m = max(l, nrow(lengthCATG))
            old = lengthCATG
            lengthCATG = matrix(0, nrow=m, ncol=ind)
            lengthCATG[1:nrow(old), ind-1] = old
            lengthCATG[l, ind] = cl
        
            m = max(s+1, nrow(scoreCATG))
            old = scoreCATG
            scoreCATG = matrix(0, nrow=m, ncol=ind)
            scoreCATG[1:nrow(old), ind-1] = old
            scoreCATG[s+1, ind] = cs

            allSeqsLengths = rep(1:nrow(lengthCATG), lengthCATG[, ind])
            allBaseScores = rep(1:nrow(scoreCATG), scoreCATG[, ind])
            dataCATG[ind, "AnalysisName"] = analysisName
            dataCATG[ind, "PicoTiterPlateID"] = picoTiterPlateID
            dataCATG[ind, "Region"] = bcm[regInd[r], 2]
            dataCATG[ind, "NumberOfReads"] = length(allSeqsLengths)
            dataCATG[ind, "NumberOfBases"] = length(allBaseScores)
            dataCATG[ind, "ReadLengthMean"] = mean(allSeqsLengths)
            dataCATG[ind, "ReadLengthStd"] = sd(allSeqsLengths)
            dataCATG[ind, "ReadLengthMedian"] = median(allSeqsLengths)
            dataCATG[ind, "ReadLengthMAD"] = mad(allSeqsLengths, constant=1)
            dataCATG[ind, "ReadScoreMean"] = mean(allBaseScores)
            dataCATG[ind, "ReadScoreStd"] = sd(allBaseScores)
            dataCATG[ind, "ReadScoreMedian"] = median(allBaseScores)
            dataCATG[ind, "ReadScoreMAD"] = mad(allBaseScores, constant=1)
            rm(allSeqsLengths)
            rm(allBaseScores)

        } else {
            stop(paste("Unknown key sequence:", bcm[keyInd[r], 2]))
        }      
    }

    names = paste(dataTCAG[,"AnalysisName"], dataTCAG[,"Region"], sep="_")
    rownames(dataTCAG) = names
    rownames(dataCATG) = names
    colnames(lengthTCAG) = names
    colnames(scoreTCAG) = names
    colnames(lengthCATG) = names
    colnames(scoreCATG) = names
    
    bcm = new("Roche454BaseCallerMetrics")
    bcm@dataTCAG = dataTCAG
    bcm@dataCATG = dataCATG
    bcm@lengthTCAG = lengthTCAG
    bcm@scoreTCAG = scoreTCAG
    bcm@lengthCATG = lengthCATG
    bcm@scoreCATG = scoreCATG

    return(bcm)
}

setMethod("read454BaseCallerMetrics", signature="character",
    .read454BaseCallerMetrics_character)
