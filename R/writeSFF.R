writeSFF <- function(sffContainer, filename) {
  
  ll <- list(flowgramFormat=flowgramFormat(sffContainer), flowChars=flowChars(sffContainer), 
             keySequence=keySequence(sffContainer), clipQualityLeft=clipQualityLeft(sffContainer), 
             clipQualityRight=clipQualityRight(sffContainer), clipAdapterLeft=clipAdapterLeft(sffContainer), 
             clipAdapterRight=clipAdapterRight(sffContainer), flowgrams=flowgrams(sffContainer), 
             flowIndexes=flowIndexes(sffContainer), reads=as.character(reads(sffContainer), use.names=TRUE), 
             qualityScores=sapply(sapply(as.character(quality(reads(sffContainer))), function(x) charToRaw(x)), 
                                  function(x) as.integer(x)-as.integer(33)))
  
  .Call("writeSFFfromR", ll, filename)
  
}
