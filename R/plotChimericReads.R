###############################################################################
## All functions below are seperated according to the standard plot, the base 
## pair plot or both The main plotting function .plotChimericReads at the end 
## calls .initializePlot, .plotReference, .annotateReference and 
## .plotReads (<-calls .plotMutations)
###############################################################################
###############################################################################

###############################################################################
################# base pair and standard plot functions #######################
###############################################################################

## .addComma helps reading (large) numbers in plots by adding commas
.addComma <- function(x){

    #paste(read.fwf(textConnection(x), commasAt, as.is = TRUE), collapse = ",")

    x = as.character(x)
    y = character()
    for(i in seq(nchar(x), 1, -3))
	if((i-3) > 0)
	    y = paste(",",substring(x,i-2,i), y, sep="")
        else
	    y = paste(substring(x,1,i), y, sep="")
    return(y)
}

## return complement strand symbol for characters ("+", or "-")
.complementStrand <- function(strand){
    if(strand == "+")
        return("-")
    if(strand == "-")
        return("+")
    stop("complementStrand: unrecognized input")
}

## sometimes the pathogenic and reciproc case are equally oriented: both show chrA <-> chrB
## this function reverses case2 with the result: pathogenic shows chrA <-> chrB, reciproc shows chrB <-> chrA
## this ensures a more intuitive visualization
.reverseC2 <- function(brpData){
    ## reverse commonBpsC2 (reference sequense, strands and local coordinates)
    refLength = length(subject(Views(commonAlignC2(brpData)[[1]])))
    commonBps = commonBpsC2(brpData)[[1]]
    tmp = commonBps$localStart
    commonBps$localStart = refLength - commonBps$localEnd + 1
    commonBps$localEnd = refLength - tmp + 1
    commonBps$strand[1] = .complementStrand(commonBps$strand[1])
    commonBps$strand[2] = .complementStrand(commonBps$strand[2])
    commonBps$refSeq = as.character(reverseComplement(DNAStringSet(commonBps$refSeq)))
    tmp = commonBps[1, ]
    commonBps[1, ] = commonBps[2, ]
    commonBps[2, ] = tmp

    ## reverse alignedReadsC2 (sequences, chromosomes, strands) 
    alnChrs = sapply(strsplit(as.character(chromosome(alignedReadsC2(brpData)[[1]])), split="/"), function(x) paste(x[2], x[1], sep="/"))
    alnStrands = sapply(strsplit(as.character(strand(alignedReadsC2(brpData)[[1]])), split="/"), function(x) paste(.complementStrand(x[1]), .complementStrand(x[2]), sep="/"))
    alnReads = reverseComplement(sread(alignedReadsC2(brpData)[[1]]))

    ## realign the reversed reads to the reversed reference sequence
    substitutionMatrix = nucleotideSubstitutionMatrix(match=1, mismatch=-1)
    alignment = pairwiseAlignment(pattern=alnReads, subject=commonBps$refSeq[1],
        type="global-local", substitutionMatrix=substitutionMatrix, gapOpening=0, gapExtension=-1)

    ## return an updated breakpoints-object(remember: the function requires, that the object holds only one breakpoint!)
    commonBpsC2(brpData) = list(commonBps)
    commonAlignC2(brpData) = list(alignment)
    alignedReadsC2(brpData) = list(AlignedRead(sread=alnReads, chromosome=factor(alnChrs), position=start(Views(alignment)), strand=factor(alnStrands)))
    return(brpData)
}

## initialize plot window
.initializePlot <- function(refStart_left, refEnd_right, plotHeight, plotBasePairs, legendColours, legendSymbols, fontScale, fontFamily){

    if(!plotBasePairs){
   	plotTol_left = 0.8		## percentage of the reference length to be added to the left plot limit
    	plotTol_right = 0.2		## percentage of the reference length to be added to the right plot limit
    }else{
    	plotTol_left = 0.3		## percentage of the reference length to be added to the left plot limit
    	plotTol_right = 0		## percentage of the reference length to be added to the right plot limit
    }

    X_min = refStart_left
    X_max = refEnd_right

    ## add some free space on the left and right limit of the plot
    xtol1 = (X_max - X_min) * plotTol_left
    xtol2 = (X_max - X_min) * plotTol_right
    X_min = X_min - xtol1
    X_max = X_max + xtol2

    plot(c(X_min, X_max), c(0, plotHeight), type="n", xaxt="n", yaxt="n", xlab="", ylab="") #, frame.plot=FALSE)

    ## add legend
    legendSize = length(legendColours)
    if(legendSize > 0){
    	legendLabels = c("deletion", "insertion", "mismatch", "breakpoint")
    	legend(x=X_min, y=plotHeight, legend=legendLabels[1:legendSize], pch=legendSymbols, col=legendColours, horiz=TRUE, cex=fontScale)
    }

    ## set character width and height to align sequences and basepairs in a basePair-plot
    if(plotBasePairs){
      characterWidth = strwidth("A", cex=fontScale, family=fontFamily)
      characterHeight = strheight("A", cex=fontScale, family=fontFamily)    
    }else{
      characterWidth = strwidth("A", cex=fontScale)
      characterHeight = strheight("A", cex=fontScale)    
    }
    return(list(characterWidth, characterHeight))

}

## annotate plot of reference sequence with chromosome, strand and position labels
.annotateReference <- function(refPos, refHeight, refStart_left, refEnd_left, refStart_right, refEnd_right, 
	chr_left, chr_right, part_left, part_right, strand_left, strand_right, labels, chrProps, fontScale, plotBasePairs){

    y = refPos + (refHeight / 2)
    ## chromosome labels
    text(x=refStart_left, y=y, label=paste("chr.", chr_left, sep=""), pos=4, cex=fontScale)
    text(x=refEnd_right, y=y, label=paste("chr.", chr_right, sep=""), pos=2, cex=fontScale)

    ## gene symbols
    text(x=refEnd_left, y=y, label=chrProps[part_left, "geneSymbol"], pos=2, cex=fontScale, font=3)
    text(x=refStart_right, y=y, label=as.character(chrProps[part_right, "geneSymbol"]), pos=4, cex=fontScale, font=3)

    ## strand labels
    ##forward
    y = refPos +  refHeight
    text(x=refStart_left, y=y, label=paste(strand_left, "5'"),pos=2, cex=fontScale)
    text(x=refEnd_right, y=y, label=paste("3'", strand_right),pos=4, cex=fontScale)
    ##reverse
    y = refPos
    if(strand_left == "+")
        text(x=refStart_left, y=y, label="- 3'",pos=2, cex=fontScale)
    else
        text(x=refStart_left, y=y, label="+ 3'",pos=2, cex=fontScale)
    if(strand_right == "+")
        text(x=refEnd_right, y=y, label="5' -",pos=4, cex=fontScale)
    else
        text(x=refEnd_right, y=y, label="5' +",pos=4, cex=fontScale)

    ## position labels (not for base pair plot)
    if(!plotBasePairs){
	offset = 0.5
    	y = refPos
    	text(x=refStart_left, y=y, pos = 1, offset=offset, label=.addComma(labels[1]), cex=fontScale)
    	text(x=refEnd_left, y=y, pos = 1, offset=offset, label=.addComma(labels[2]), cex=fontScale)
    	y = refPos + refHeight
    	text(x=refStart_right, pos = 3, offset=offset, y=y, label=.addComma(labels[3]), cex=fontScale)
    	text(x=refEnd_right, pos = 3, offset=offset, y=y, label=.addComma(labels[4]), cex=fontScale)		
    }
}

## .getMutations reads indels and mismatches from the alignment and modifies the coordinates for correct plotting
.getMutations <- function(alignment, from, to){

    alnStart = start(Views(alignment)) - 1

    ## deletions, insertions and mismatches from alignment
    dels = deletion(alignment)[[1]]
    ins = insertion(alignment)[[1]]
    mms = mismatchTable(alignment)

    ## adjust coordinates of mismatches and deletions according to insertions (and deletions) removed prior to them
    if(nrow(mms) > 0)
        for(m in 1:nrow(mms)){
	    mms[m, c("PatternStart", "PatternEnd")] = 
		mms[m, c("PatternStart", "PatternEnd")] - sum(width(restrict(ins, 1, mms$PatternStart[m]))) +
		    sum(width(restrict(dels, 1, mms$PatternStart[m])))
	}
    if(length(dels) > 0)
        for(d in 1:length(dels))
          dels[d] = shift(dels[d], -sum(width(restrict(ins, 1, start(dels[d])))))

    ## adjust coordinates of dels/ins/mismatches according to the beginning of alignment on the reference sequence
    ## take care of dels ("-") inserted prior to a deletion
    if(length(dels) > 0)
      dels = shift(dels, seq(0, length(dels) - 1) + alnStart)
    dels = as.data.frame(dels)
    ins = shift(ins, alnStart)
    ins = as.data.frame(ins)
    mms[ , c("PatternStart", "PatternEnd")] = mms[ , c("PatternStart", "PatternEnd")] + alnStart
    ## only draw dels/ins/mismatches within a given range (from/to; +/- maxBasePairs around the breakpoint)
    dels = subset(dels, start >= from & end <= to)
    ins = subset(ins, start >= from & end <= to)
    mms = subset(mms, PatternStart >= from & PatternEnd <= to)

    if(nrow(dels) > 0)
    	dels = unlist(c(apply(dels, 1, function(x) x["start"]:x["end"]))) - (from - 1)
    else
	dels=c()
    if(nrow(ins) > 0)
    	ins = ins$start - (from - 1)
    else
	ins=c()
    if(nrow(mms) > 0)
    	mms = mms$PatternStart - (from - 1)
    else
	mms=c()

    return(list(dels, ins, mms))

}

###############################################################################
################# standard plot functions #####################################
###############################################################################

## draw reference sequence in the center of the plot
.plotReference <- function(refHeight, refPos, refStart_left, refStart_right, refEnd_left, refEnd_right, 
	part_left, part_right, chrProps, colours){

    y = c(refPos, refPos + refHeight, refPos + refHeight, refPos)

    ## left part
    x = c(refStart_left, refStart_left, refEnd_left, refEnd_left)
    polygon(x, y, col=chrProps[part_left, "colourRef"])
    ## right part
    x = c(refStart_right, refStart_right, refEnd_right, refEnd_right)
    polygon(x, y, col=chrProps[part_right, "colourRef"])

    ## gap between both parts
    x = c(refEnd_left, refEnd_left, refStart_right, refStart_right)
    polygon(x, y, density=c(20,30), angle=c(-45,45), col=colours[10], border=colours[10])
}

## plots deletions, insertions and mismatches for standard plot
.plotMutations = function (dels, ins, mms, y, colours){

    points(x=dels, y=rep(y,length(dels)), pch="|", col=colours[7])	
    points(x=ins, y=rep(y,length(ins)), pch="|", col=colours[8])	
    points(x=mms, y=rep(y,length(mms)), pch="|", col=colours[9])	

}

## draw reads above and below the reference sequence
.plotReads <- function(alignments, names, directions, numReads, refPos, refHeight, readDist, localBrp_left, localBrp_right,
	refEnd_right, refStart_left, plotMut, colours, fontScale, characterHeight){

    yDist_above = 1
    yDist_below = 1

    for(i in 1:numReads){

   	## x-positions
	start = start(Views(alignments[i]))
	end = end(Views(alignments[i]))

    	## y-position (depends on direction: distinguish between two directions: 3'-5' = below reference, 5'-3' = above reference)
	## forward (5'-3')
       	if(directions[i] == "forward"){
    	    y = refPos + refHeight + 3*characterHeight + (yDist_above * readDist)
    	    yDist_above = yDist_above + 1	
	## reverse (3'-5')
	}else{
    	    y = refPos - 3*characterHeight - (yDist_below * readDist)
    	    yDist_below = yDist_below + 1	
	}

	## plot read and mutations
    	lines(x=c(start, end), y=c(y, y), col=colours[6])
	points(x=start, y=y, pch=20, col=colours[6])
	points(x=end, y=y, pch=20, col=colours[6])
	points(x=c(localBrp_left, localBrp_right - 1), y=c(y, y), pch="|", col=colours[10])

    	if(plotMut){
    	    muts = .getMutations(alignments[i], refStart_left, refEnd_right)
       	    .plotMutations(muts[[1]], muts[[2]], muts[[3]], y, colours)
	}

	## label each read with its name
	labelOffset = 2
      	    text(x=refStart_left, y=y, label=names[i], pos=2, offset = labelOffset, cex=fontScale)
	
    }	
}

###############################################################################
################# base pair plot functions ####################################
###############################################################################

## draw reference sequence in the center of the plot
.plotReferenceBasePairs <- function(numReadsAbove, numReadsBelow, readDist, refHeight, refPos, refSeq, localBrp_left, localBrp_right,
	refStart_left, refStart_right, refEnd_left, refEnd_right, part_left, part_right, from, to, chrProps, colours, labels, 
	fontScale, fontFamily, characterWidth, characterHeight, maxBasePairs_left, maxBasePairs_right){

    x_left = c(refStart_left, refStart_left, refEnd_left,  refEnd_left)
    x_right = c(refStart_right, refStart_right, refEnd_right, refEnd_right)
    x_gap = c(refEnd_left, refEnd_left, refStart_right, refStart_right)
    yReads = c(refPos - 2*characterHeight - readDist * numReadsBelow, refPos + refHeight + 3*characterHeight + readDist * numReadsAbove,
	refPos + refHeight + 3*characterHeight + readDist * numReadsAbove, refPos - 2*characterHeight - readDist * numReadsBelow)
    yRef = c(refPos - characterHeight,refPos + refHeight + characterHeight,refPos + refHeight + characterHeight,refPos - characterHeight)

    ## draw coloured frames for reads (left and right)
    polygon(x_left, yReads, col=chrProps[part_left, "colourReads"], border="black")
    polygon(x_right, yReads, col=chrProps[part_right, "colourReads"], border="black")

    ## draw coloured frames for reference sequence (left and right)
    polygon(x_left, yRef, col=chrProps[part_left, "colourRef"], border="black")
    polygon(x_right, yRef, col=chrProps[part_right, "colourRef"], border="black")

    ## draw coloured frame for the gap
    polygon(x_gap, yReads, col=colours[5], border=colours[5])

    rect(xleft=refStart_left, ybottom=refPos - 2*characterHeight - readDist * numReadsBelow,
   	xright=refEnd_right, ytop=refPos + refHeight + 3*characterHeight + readDist * numReadsAbove, col=NA, border="black")

    ## draw basepairs of the reference sequence; compute reverse/reverseComplement for 3'-5' direction
    text(x=refStart_left, y=refPos + refHeight, label=substr(refSeq, from, to), pos=4, offset=0, cex=fontScale, family=fontFamily, col=colours[6])
    text(x=refStart_left, y=refPos, label=substr(complement(DNAString(refSeq)), from, to), pos=4, offset=0, cex=fontScale, family=fontFamily, col=colours[6])
}

## plots deletions and insertions for basePair-plot
.plotMutationsBasePairPlot <- function(x, y, dels, ins, mms, fontScale, fontFamily, characterWidth, characterHeight, colours){

    ## mark mutations as triangles
    mutHeight = characterHeight * 0.6
    mutWidth = characterWidth * 0.6
    y = y + characterHeight

    ## print deletions as (red) triangle 
    if(length(dels) > 0){
    	dels = dels * characterWidth
    	xd = x - 0.5 * characterWidth  ## draw in the center of a symbol
    	sapply(dels, function(d) 
	    polygon(x=c(xd + d, xd + d - mutWidth, xd + d + mutWidth), y=c(y, y + mutHeight, 
	    y + mutHeight), col=colours[7], border=colours[7]) 
	)
    }
    ## print insertions as (green) triangle
    if(length(ins) > 0){
    	ins = ins * characterWidth
	xi = x - characterWidth
    	sapply(ins, function(i) 
	    polygon(x=c(xi + i, xi + i - mutWidth, xi + i + mutWidth), y=c(y, y + mutHeight, 
	    y + mutHeight), col=colours[8], border=colours[8]) 
	)
    }
    ## print matches as (black) triangle
    if(length(mms) > 0){
    	mms = mms * characterWidth
    	xm = x - 0.5 * characterWidth  ## draw in the center of a symbol
    	sapply(mms, function(m) 
	    polygon(x=c(xm + m, xm + m - mutWidth, xm + m + mutWidth), y=c(y, y + mutHeight, 
	    y + mutHeight), col=colours[9], border=colours[9]) 
	)
    }
}

## draw reads above and below the reference sequence
.plotReadsBasePairs <- function(alignments, names, directions, numReads, refPos, refHeight, readDist, refStart_left, localBrp_left, localBrp_right, 
	from, to, plotMut, colours, fontScale, fontFamily, characterWidth, characterHeight){

    yDist_above = 1
    yDist_below = 1

    for(i in 1:numReads){

	seq = toString(aligned(alignments)[i])

   	## x-position
        x = refStart_left

    	## y-position (depends on direction: distinguish between two directions: 3'-5' = below reference, 5'-3' = above reference)
	## forward (5'-3')
       	if(directions[i] == "forward"){
    	    y = refPos + refHeight + characterHeight + (yDist_above * readDist)
    	    yDist_above = yDist_above + 1	
	## reverse (3'-5')
	}else{
    	    y = refPos - characterHeight - (yDist_below * readDist)
    	    yDist_below = yDist_below + 1	
	    seq = toString(complement(DNAString(seq)))
	}

	## only plot basepairs within a given range (from, to; +/- maxBasePairs around the breakpoint)
	seq = substr(seq, from, to)

	## plot read and mutations
    	text(x=x, y=y, label= seq, pos=4, offset=0, cex=fontScale, family=fontFamily, col=colours[6])
    	if(plotMut){
    	    muts = .getMutations(alignments[i], from, to)
	    .plotMutationsBasePairPlot(x, y, muts[[1]], muts[[2]], muts[[3]], fontScale, fontFamily, characterWidth, characterHeight, colours)	    
	}

	## label each reads with its name
	labelOffset = 2
       	text(x=x, y=y, label=names[i], pos=2, offset = labelOffset, cex=fontScale)
	
    }	
}

##############################################################################
################# Main function ##############################################
##############################################################################

.plotChimericReads <- function(brpData, geneSymbols=FALSE, plotMut=TRUE, plotBasePairs=FALSE, maxBasePairs=50, legend=FALSE, title="", 
    col=c("red", "green", "black", "orange")){

    ## catch some bad input
    if(length(brpData) > 1)
	stop("invalid argument: please pass just one breakpoint cluster for plotting (try subsetting your breakpoint data)")
    if(!is.logical(plotMut))
	stop("invalid argument: please specify a logical value for plotMut")    
    if(!is.logical(plotBasePairs))
	stop("invalid argument: please specify a logical value for plotBasePairs")
    if(!is.logical(geneSymbols))
	if(!is.vector(geneSymbols, mode="character") | (length(geneSymbols) != 2))
	    stop("invalid argument: please specify a logical value or a vector of two strings for geneSymbols")
    if(!is.numeric(maxBasePairs) | (maxBasePairs < 1))
	stop("invalid argument: please specify value > 0 for maxBasePairs")
    if(maxBasePairs > 1000)
	warning(paste("sure you want to plot up to", 2*maxBasePairs, "base pairs of each read?"))
    if(!is.character(title))
	stop("invalid argument: title has to be of type character")
    if(!is.vector(col, mode="character") | (length(col) != 4))
	stop("invalid argument: please specify four valid colours for deletions, insertions, mismatches and breakpoints")

    ## some (heuristic) constant variables:
    ## scale and font family (font family only for basepairs in basePair-plot)
    fontScale = 0.7
    fontFamily = "mono"
    ## distance between reads, height of the reference sequence 
    if(!plotBasePairs){
      	readDistInit = 4
       	refHeight = 5
    }else{
    	readDistInit = 0.1	## distance of reads for initialization of plot window (later depending on character height, but not known yet)
    	refHeight = 0.1
    }
    ## maximum distance to a breakpoint say it belongs to another one (only needed to distinguish between both parts of an inversion)
    maxInvDist = 1000

    ## define plot colours
    ## in this order: reference background (left), reads background (left), reference background (right), reads background (right),
    ## breakpoints (standard plot only), gap (base pair plot only), reference/reads
    ## colors for deletions, insertions, mismatches and breqakpoints given by function argument (default: red, green, black, orange)
    colours = c("lightsteelblue", "lightsteelblue1", "lightyellow2", "ivory", "lightgrey", "black", col)


    ## initialize biomaRt
    if(class(geneSymbols) == "logical")
	if(geneSymbols)
    	    ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

    ## see if data requires two plots (case2 contains breakpoint and alignment data)
    if(nrow(commonBpsC2(brpData)[[1]]) == 0)
      numPlots = 1
    else
      numPlots = 2

    ## split plot window into two parts if necessary (horizontally for basepair-plot vertically for standard plot)
    if(numPlots == 2)
    	if(plotBasePairs)
    	    par(mfrow=c(2, 1), oma=c(0, 0, 2, 0), mar = c(0.5, 0.5, 0.5, 0.5))
    	else
            ## for horizontal alignment    
    	    par(mfrow=c(1,2), oma=c(0, 0, 2, 0), mar = c(0.5, 0.5, 0.5, 0.5))
#            ## for vertical alignment    
#    	    par(mfrow=c(2,1), oma=c(0, 0, 2, 0), mar = c(0.5, 0.5, 0.5, 0.5))
    else
    	par(mfrow=c(1,1), oma=c(0, 0, 2, 0), mar = c(0.5, 0.5, 0.5, 0.5))

    ## if the pathogenic and reciproc case are equally oriented (both show chrA <-> chrB), reverse case2 (reciproc case)
    ## in case of an inversion (chrA == chrB), compare the first two breakpoints instead of the first two chromosomes
    if(numPlots == 2){
    	if(commonBpsC1(brpData)[[1]]$chr[1] != commonBpsC1(brpData)[[1]]$chr[2]){
	    if(commonBpsC1(brpData)[[1]]$chr[1] == commonBpsC2(brpData)[[1]]$chr[1])
            	brpData = .reverseC2(brpData)
   	}else{
	    if(abs(commonBpsC1(brpData)[[1]]$breakpoint[1] - commonBpsC2(brpData)[[1]]$breakpoint[1]) < 100)
            	brpData = .reverseC2(brpData)
    	}
    }

    ## retrieve read names and directions (forward: 5'-3', reverse: 3'-5')
    ## according its the direction a read will be plotted below or above the reference sequence
    ##-> count the reads below and above the reference sequence to correctly initialize the plot window
    seqNames_all = list(length=numPlots)
    directions_all = list(length=numPlots)
    numReadsAbove_all = vector(mode="numeric", length=numPlots)
    numReadsBelow_all = vector(mode="numeric", length=numPlots)

    refStrands = paste(commonBpsC1(brpData)[[1]]$strand, collapse="/") ## note: strand information is formatted e.g. as "+/-" in class AlignedRead
    seqNames_all[[1]] = names(sread(alignedReadsC1(brpData)[[1]]))
    directions_all[[1]] = sapply(strand(alignedReadsC1(brpData)[[1]]), function(x) if(x != refStrands) return("reverse") else return("forward")) 
    numReadsAbove_all[1] = sum(directions_all[[1]] == "forward")
    numReadsBelow_all[1] = sum(directions_all[[1]] == "reverse")
    if(numPlots == 2){
    	refStrands = paste(commonBpsC2(brpData)[[1]]$strand, collapse="/") ## note: strand information is formatted e.g. as "+/-" in class AlignedRead
    	seqNames_all[[2]] = names(sread(alignedReadsC2(brpData)[[1]]))
    	directions_all[[2]] = sapply(strand(alignedReadsC2(brpData)[[1]]), function(x) if(x != refStrands) return("reverse") else return("forward")) 
    	numReadsAbove_all[2] = sum(directions_all[[2]] == "forward")
    	numReadsBelow_all[2] = sum(directions_all[[2]] == "reverse")
    }
    numReadsAbove_max = max(numReadsAbove_all)
    numReadsBelow_max = max(numReadsBelow_all)

    ## plot reference sequence and reads
    for(d in 1:numPlots){

	## get information from breakpoint data (list: 1.) breakpoints, 2.) alignments, 3.) reads)
    	if(d == 1){
	    breakpoints = commonBpsC1(brpData)[[1]]
 	    alignment = commonAlignC1(brpData)[[1]]
	}else{
	    breakpoints = commonBpsC2(brpData)[[1]]
 	    alignment = commonAlignC2(brpData)[[1]]
	}
    	numReads = length(alignment)

	## initialize all information about the reads (chromosomes and strands etc.)
	chr_left = part_left = breakpoints$chr[1]
    	chr_right = part_right = breakpoints$chr[2]
    	strand_left = breakpoints$strand[1]
    	strand_right = breakpoints$strand[2]
	geneSymbol_left = ""
	geneSymbol_right = ""
	seqNames = seqNames_all[[d]]
	directions = directions_all[[d]]
	numReadsAbove = numReadsAbove_all[[d]]
	numReadsBelow = numReadsBelow_all[[d]]

	## initialize all information about the reference (breakpoint etc.)
	refSeq = breakpoints$refSeq[1]
	brp_left = breakpoints$breakpoint[1]
	brp_right = breakpoints$breakpoint[2]
	localBrp_left = breakpoints$localEnd[1]
	localBrp_right = breakpoints$localStart[2]
	refLength_left = localBrp_left
	refLength_right = breakpoints$localEnd[2] - localBrp_right + 1

    	## set position labels (from original alignment data)
	if(strand_left == "+")
    	    label_start = brp_left - refLength_left
	else
    	    label_start = brp_left + refLength_left
	if(strand_right == "+")
    	    label_end = brp_right + refLength_right
	else
    	   label_end = brp_right - refLength_right
    	labels = c(label_start, brp_left, brp_right, label_end)

    	## load gene symbols from biomaRt
	if(class(geneSymbols) != "logical"){
	    geneSymbol_left = geneSymbols[1]
	    geneSymbol_right = geneSymbols[2]
	}else{
	    if(geneSymbols){
    	    	geneSymbol_left = getBM(attributes=c("hgnc_symbol"), 
        	    filter = c("chromosome_name", "start", "end"), 
		    values = list(chr_left, labels[2], labels[2]), 
		    mart = ensembl)
    	    	geneSymbol_right = getBM(attributes=c("hgnc_symbol"), 
        	    filter = c("chromosome_name", "start", "end"), 
		    values = list(chr_right, labels[3], labels[3]), 
		    mart=ensembl)
		## paste gene symbols if more than one reported
		geneSymbol_left = paste(geneSymbol_left$hgnc_symbol, collapse=",")
		geneSymbol_right = paste(geneSymbol_right$hgnc_symbol, collapse=",")
	    }
	}
	## assign colours and geneSymbols to the two parts (chromosomes) of the chimeric reads (only once for first plot)
	## in case of inversion (equal chromosome): assign a name (part1/2) to each part according to location of the breakpoint
        if(d == 1){
	    if(chr_left == chr_right){
		brp_part1 = brp_left
		part_left = "part1"
		part_right= "part2"
	    }
            chrProps = data.frame(
	        colourRef=c(colours[1], colours[3]), 
	        colourReads=c(colours[2], colours[4]),
	        geneSymbol=c(geneSymbol_left, geneSymbol_right),
	        row.names=c(part_left, part_right),
	        stringsAsFactors=FALSE
	    )
	}else{
	    if(chr_left == chr_right)
	        if(abs(brp_left - brp_part1) < maxInvDist){
	            part_left = "part1"
	    	    part_right= "part2"
		}else{
		    part_left = "part2"
		    part_right= "part1"
		}
   	}

	## initialize coordinates of the reference sequence
	refStart_left = 1
	refEnd_left = localBrp_left
	refStart_right = localBrp_right - 1
	refEnd_right = nchar(refSeq)

	## for base pair plots: only plot a subset of the base pairs around th breakpoint
    	maxBasePairs_left = min(maxBasePairs, refLength_left)
    	maxBasePairs_right = min(maxBasePairs, refLength_right)

    	from = localBrp_left - maxBasePairs_left + 1
    	to = localBrp_right + maxBasePairs_right - 1

	## initialize plot window
	if(!plotBasePairs){
    	    plotHeight = ((numReadsAbove_max + numReadsBelow_max + 2) * readDistInit) + refHeight
	    refPos = ((numReadsBelow_max + 1) * readDistInit) - (refHeight / 2)
	    ## colours for the legend: colours for dels, ins, mismatches and brp (legend only for first plot)
	    if(!legend | d == 2) 
		legendColours = c()
	    else 
		legendColours = colours[7:10]
	    legendSymbols = rep("|", 4)
	    init = .initializePlot(refStart_left, refEnd_right, plotHeight, plotBasePairs, legendColours, legendSymbols, fontScale, fontFamily)
	}
	else{
    	    plotHeight = ((numReadsAbove + numReadsBelow + 2) * readDistInit) + refHeight
	    refPos = ((numReadsBelow + 1) * readDistInit) - (refHeight / 2)
	    if(!legend | d == 2) 
		legendColours = c()
	    else 
		legendColours = colours[7:9]
	    legendSymbols = rep(6, 3)	## triangle
	    init = .initializePlot(refStart_left, refStart_right + maxBasePairs_right, plotHeight, plotBasePairs, legendColours, 
		legendSymbols, fontScale, fontFamily)
	}
    	characterWidth = init[[1]]
	characterHeight = init[[2]]
	readDist = characterHeight * 3

	## when plotting base pairs, adjust reference coordinates according to the character widths and plotting region
	if(plotBasePairs){
	    gap = refStart_right - refEnd_left
	    refEnd_left = refStart_left + (characterWidth * maxBasePairs_left)	## ends after symbol at poisiton maxBasePairs
	    refStart_right = refEnd_left + (characterWidth * gap)			## starts at symbol at poisiton maxBasePairs + 1
	    refEnd_right = refStart_right + (characterWidth * maxBasePairs_right)  	## ends after last symbol at poisiton 2*maxBasePairs
	}

	## draw reference
    	if(!plotBasePairs)
	    .plotReference(refHeight, refPos, refStart_left, refStart_right, refEnd_left, refEnd_right, part_left, part_right, chrProps, colours)
	else
	    .plotReferenceBasePairs(numReadsAbove, numReadsBelow, readDist, refHeight, refPos, refSeq, localBrp_left, localBrp_right,
		refStart_left, refStart_right, refEnd_left, refEnd_right, part_left, part_right, from, to, chrProps, colours, labels, 
		fontScale, fontFamily, characterWidth, characterHeight, maxBasePairs_left, maxBasePairs_right)
	.annotateReference(refPos, refHeight, refStart_left, refEnd_left, refStart_right, refEnd_right, 
	    chr_left, chr_right, part_left, part_right, strand_left, strand_right, labels, chrProps, fontScale, plotBasePairs)

    	## draw reads
    	if(!plotBasePairs)
	    .plotReads(alignment, seqNames, directions, numReads, refPos, refHeight, readDist, localBrp_left, localBrp_right,
	        refEnd_right, refStart_left, plotMut, colours, fontScale, characterHeight)
	else
	    .plotReadsBasePairs(alignment, seqNames, directions, numReads, refPos, refHeight, readDist, refStart_left, 
		localBrp_left, localBrp_right, from, to, plotMut, colours, fontScale, fontFamily, characterWidth, characterHeight)

    }

    ## Finally, the plot receives its title (if any)
    mtext(title, side=3, outer=TRUE) 
}

setMethod("plotChimericReads",
    signature(brpData="Breakpoints"),
    .plotChimericReads
)
