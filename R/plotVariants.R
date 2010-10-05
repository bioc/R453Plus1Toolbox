#########################################################################################
## plotVariants creates a plot of all variant positions in a given transcript.
## input:		
##  variant data as returned by annotateVariants,
##  bmInfo ensembl transcript information,
##  legend TRUE/FALSE (optional),
##  regions to be shaded within the plot (optional),
##  region label for the legend (optional),
##  col four colours for missense, silent, nonsense mutations and deletions
##  title (optional),
## output: nothing but a beautiful plot
#########################################################################################

.plotVariants <- function(varData, bmInfo, legend, regions, regionsGiven, regLabel, col, title){


  message("Creating variant plot for transcript ", bmInfo$ensembl_transcript_id[1], " ... ", appendLF = FALSE)
  transcript = bmInfo$ensembl_transcript_id[1]
  exon_end = bmInfo$cds_length[1]	
  
  ## catch more bad input:
  if(exon_end %% 3 == 0)
    exon_end = exon_end / 3
  else
    stop("problem with transcript length")
  
  ## regions (if any) have to consists of start and end points within a suitable range (1-exon_end)
  if(regionsGiven)
    if(any(regions > exon_end))
      stop("given regions out of range")

  ## draw basic transcript structure:
  par(mar = c(0.5, 0.5, 1, 0.5))
  plot(c(1,exon_end), c(0,50), type="n", xaxt="n", yaxt="n", xlab="",
       ylab="", frame.plot=FALSE)
  title(title)
  
  ## mark given regions
  if(regionsGiven){
    for(i in seq(1,length(regions),2)){
      x = c(regions[i], regions[i], regions[i+1], regions[i+1])
      y = c(25,40,40,25)
      polygon(x, y, col="gray", border="gray", density=c(20,30),
              angle=c(-45,45)) 
    }
  }
  
  ## draw exon regions
  for(i in 1:nrow(bmInfo)){
    start = bmInfo[i,"cds_start"] / 3
    end = bmInfo[i,"cds_end"] / 3
    x = c(start, start, end, end)
    y = c(25 ,40, 40, 25)
    polygon(x,y)
    text((start + end) / 2, 28, bmInfo[i, "rank"])
  }
  axis(1, c(1, seq(100, exon_end,100), exon_end), pos=15)
  text(x=exon_end/2, y=10, "position (amino acid)", pos=1)
  
  ## draw variant positions
  ## calculate the drawing width according to the codon length
  characterWidth = strwidth("A", family="mono")
  characterHeight = strheight("A", family="mono")
  width = characterWidth
  height = 3
  yMisMut = 40
  yNonMut = 40 - height
  yDel = 25 - height	
  ## for the legend: remember if a mutation appears in the plot
  nonMutPresent = FALSE
  misMutPresent = FALSE
  silentMutPresent = FALSE
  delPresent = FALSE
  numNonMut = 0
  numMisMut = 0
  numSilentMut = 0
  numDel = 0

  varData = annotatedVariants(varData)
  varNames = names(varData)
  ## test every variant if it lies within the transcript
  for(vn in varNames){
    
    v = varData[[vn]]$exons
    ## see if exon data available
    if(nrow(v) > 0){
      
      ## see if variant lies on the given transcript
      v = subset(v, v$ensembl_transcript_id == transcript)
      if(nrow(v) > 0){

        for(i in 1:nrow(v)){
          
          ## draw markers for variants that lie in a coding area
          if(v[i,"coding"] == TRUE){
            if(length(grep("-", v[i, "codonMut"])) == 0){
              pos = v[i, "numCodonStart"]
              
              if(v[i,"AminoRef"] != "*" & v[i,"AminoMut"] != "*"){
                
                ## white marker for missense point mutation
                if(v[i,"AminoRef"] !=  v[i,"AminoMut"]){
                  .plotMisMut(pos, yMisMut, width, height, col[1])
                  misMutPresent = TRUE
                  numMisMut = numMisMut + 1
                }else{
                  ## grey marker for silent point mutation
                  .plotSilentMut(pos, yMisMut, width, height, col[2])
                  silentMutPresent = TRUE
                  numSilentMut = numSilentMut + 1
                }
                
                ## black marker for nonsense point mut
              } else {
                if(v[i,"AminoRef"] != "*" & v[i,"AminoMut"] == "*"){
                  .plotNonMut(pos, yNonMut, width, height, col[3])
                  nonMutPresent = TRUE
                  numNonMut = numNonMut + 1
                }else
                warning("variant ", vn, " is of unknown type")
                
              }
              
              ## white diamond marker for deletion
            } else { 
              pos = (v[i, "numCodonStart"] +
                     v[i, "numCodonEnd"]) / 2
              .plotDel(pos, yDel, width, height, col[4])
              delPresent = TRUE
              numDel = numDel + 1
            }
          }
        }
      }
    }
  }
  
  ## draw legend
  if(legend){
    
    nonMutLabel = "nonsense mutation"
    misMutLabel = "missense mutation"
    silentMutLabel = "silent mutation"
    delLabel = "deletion"
    
    nonMutLength = nonMutPresent * (1 + 4*width + nchar(nonMutLabel)*characterWidth) ## zero, if no such mutation present
    misMutLength = misMutPresent * (1 + 4*width + nchar(misMutLabel)*characterWidth)
    silentMutLength = silentMutPresent * (1 + 4*width + nchar(silentMutLabel)*characterWidth)
    delLength = delPresent * (1 + 4*width + nchar(delLabel)*characterWidth)
    x = c(1, 1 + nonMutLength, 1 + nonMutLength + misMutLength,  1 + nonMutLength + misMutLength + silentMutLength, 1 + nonMutLength + misMutLength + silentMutLength +delLength)
    x[c(!nonMutPresent, !misMutPresent, !silentMutPresent, !delPresent, !regionsGiven)] = NA		## omit labels that do not appear in the plot
    y = 0
    
    
    ## nonsense mutation label
    .plotNonMut(x[1], y, width, height, col[3])
    text(x[1] + 2*width, y + (characterHeight/2), nonMutLabel, pos=4)
    
    ## missense mutation label
    .plotMisMut(x[2], y, width, height, col[1])
    text(x[2] + 2*width, y + (characterHeight/2), misMutLabel, pos=4)
    
    ## silent mutation label
    .plotSilentMut(x[3], y, width, height, col[2])
    text(x[3] + 2*width, y + (characterHeight/2), silentMutLabel, pos=4)
    
    ## del mutation label
    .plotDel(x[4], y, width, height, col[4])
    text(x[4] + 2*width, y + (characterHeight/2), delLabel, pos=4)

    ## region label
    xReg = c(x[5] - width, x[5] -width, x[5] + width, x[5] + width)
    yReg = c(y, y + height, y + height, y)
    polygon(xReg, yReg, col="gray", border="gray", density=c(20,30), angle=c(-45,45)) 
    text(x[5] + 2*width, y + (characterHeight/2), regLabel, pos=4)
  }

  message("done\n")
  message("Missense mutations: ", numMisMut)
  message("Nonsense mutations: ", numNonMut)
  message("Silent mutations: ", numSilentMut)
  message("Deletions: ", numDel)
}

## additional functions to draw mutations
## defaults:
## nonsense mutations as black triangles
## missense mutations as grey triangles
## silent mutations as white triangles
## deletions as white diamonds
.plotNonMut <- function(x, y, width, height, col){
  x = c(x, x - width, x + width)	
  y = c(y + height, y, y)
  polygon(x, y, col=col)
}
.plotMisMut <- function(x, y, width, height, col){
  x = c(x, x - width, x + width)
  y = c(y, y + height, y + height)
  polygon(x, y, col=col)
}
.plotSilentMut <- function(x, y, width, height, col){
  x = c(x, x - width, x + width)
  y = c(y, y + height, y + height)
  polygon(x, y, col=col)
}
.plotDel <- function(x, y, width, height, col){
  x = c(x, x - width, x, x + width)
  y = c(y, y + (height/2), y + height, y + (height/2))
  polygon(x, y, col=col)
}


.checkInput <- function(title, legend, regions, regionsGiven, regLabel, col){
  
  ## catch some bad input:
  if(!is.character(title))
    stop("invalid argument: title has to be of type character")
  if(!is.character(regLabel))
    stop("invalid argument: regLabel has to be of type character")
  if(!is.logical(legend))
    stop("invalid argument: please specify a logical value (TRUE/FALSE) for legend")
  if(!is.vector(col, mode="character") | (length(col) != 4))
    stop("invalid argument: please specify four valid colours for missense, silent, nonsense mutations and deletions") 
  
  ## regions (if any) have to consists of start and end points within a suitable range (1-exon_end)
  regionsGiven = length(regions) > 0
  if(regionsGiven)
    if(!is.vector(regions, mode="numeric") || length(regions) %% 2 != 0 || any(regions < 1))
      stop("invalid argument: given regions in wrong format; please specify start and end position for every region within a suitable range")
  
}


setMethod("plotVariants",
          signature=c(varData="AnnotatedVariants", transcript="character"),
          
          function(varData, transcript, legend=TRUE, regions=c(), regLabel="(region label)", col=c("grey","white","black","white"), title=""){
            
            message("Retrieving coding regions from Ensembl for transcript ", transcript, ":")
            regionsGiven = length(regions) > 0
            .checkInput(title, legend, regions, regionsGiven, regLabel, col)
            ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")
            bmAttributes=c(
              "ensembl_transcript_id",
              "rank",
              "cds_start",
              "cds_end",
              "cds_length"
              )
            bmInfo = getBM(attributes=bmAttributes, filters="ensembl_transcript_id",
              values=transcript, mart=ensembl)
            
            ## only select coding exons
            bmInfo = subset(bmInfo, !(is.na(bmInfo$cds_start) & is.na(bmInfo$cds_end)))
            .plotVariants(varData, bmInfo, legend, regions, regionsGiven, regLabel, col, title)
            return(bmInfo)
            
          })


setMethod("plotVariants",
          signature=c(varData="AnnotatedVariants", transcript="data.frame"),
          
          
          function(varData, transcript, legend=TRUE, regions=c(), regLabel="(region label)", col=c("grey","white","black","white"), title=""){
            
            if(!all(c("ensembl_transcript_id","rank","cds_start","cds_end","cds_length") %in% colnames(transcript)))
              stop("invalid input: transcript information requires at least the following colummns: ensembl_transcript_id, rank, cds_start, cds_end, cds_length (try calling plotVariants using an ensembl transcript-id)")
            regionsGiven = length(regions) > 0
            .checkInput(title, legend, regions, regionsGiven, regLabel, col)
            
            .plotVariants(varData, transcript, legend, regions, regionsGiven, regLabel, col, title)
            
          })
