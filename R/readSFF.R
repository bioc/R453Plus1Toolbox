#dyn.load("readSFF.so")

readSFF = function(files) {

  # check or create the list of files
  if (file.exists(files[1])) {
    tmp = file.info(files[1])
    if (tmp$isdir) { # file is directory
      directory <- files[1]
      directory <- sub("[\\/]$", "", directory)
      sfffiles <- try (dir(directory, full.names=TRUE))
      if (inherits(sfffiles, "try-error") || length(sfffiles) == 0) stop("Can't read directory: ", directory)
      sfffiles <- sfffiles[tolower(substring(sfffiles, nchar(sfffiles)-3)) == ".sff"]
    } else { # read invididual files	
      sfffiles <- files
    }
  } else {
    stop("Please provide a valid file or directory name.")
  }	

  containerlist = list()
  for(file in sfffiles) {
    if(!file.exists(file)) {
      cat("Error reading file: ", file ," - skipping.\n", sep="")
      next
    }

    cat("Reading file ", basename(file), " ...", sep="")
    result = .Call("readSFFfromR", file, PACKAGE="R453Plus1Toolbox")
  
    sffcontainer = new("SFFContainer")
  
    filename(sffcontainer) = basename(file)
    flowgramFormat(sffcontainer) = result[["flowgramFormat"]]
    flowChars(sffcontainer) = result[["flowChars"]]
    keySequence(sffcontainer) = result[["keySequence"]]
  
    clipQualityLeft(sffcontainer) = result[["clipQualityLeft"]]
    clipQualityRight(sffcontainer) = result[["clipQualityRight"]]
    clipAdapterLeft(sffcontainer) = result[["clipAdapterLeft"]]
    clipAdapterRight(sffcontainer) = result[["clipAdapterRight"]]
    flowgrams(sffcontainer) = result[["flowgrams"]]
    flowIndexes(sffcontainer) = result[["flowIndexes"]]
    qualityScores = PhredQuality(sapply(result[["qualityScores"]], function(x) rawToChar(as.raw(x + 33)))) 
    reads(sffcontainer) = QualityScaledDNAStringSet(DNAStringSet(result[["reads"]]), qualityScores)
  
    containerlist[[basename(file)]] = sffcontainer 
    cat(" done! \n", sep="")
  }

  if(length(containerlist) > 1) {
    return(containerlist)
  } else {
    return(containerlist[[1]])
  }
  
}
