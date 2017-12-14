
.plotVariants <- function(data, gene, transcript, bmExons, bmGene, bmReturn, regions, mutationInfo, groupBy, horiz, cex, title, legend){

  ## When sorting by position, take care that different mutations at the same position will not be in the same group
  ## Therefore, combine the position with the mutation type
  if(groupBy == "pos"){
    groupBy = "groups"
    data$groups = paste(data$pos, data$mutation, sep="_")
  }
  data = data[order(data[, groupBy]), ]
  
  ## parameter horiz stands for horizontal alignment of a group of mutations a single position
  ## if horiz is false, these mutations will be plotted vertically
  if(horiz == FALSE){
    ## list of mutations for every group
    mutGroup = !duplicated(data[, groupBy])
    muts = list()
    for(i in 1:sum(mutGroup)){
      g = data[mutGroup, groupBy][i]
      muts[[i]] = as.character(data[data[, groupBy]==g, "mutation"])
    }
    numMuts = length(unique(data[ , groupBy]))
  }else{
    ## if horiz is true every mutation will be plotted horizontally
    mutGroup = 1:nrow(data)
    muts = as.list(data$mutation)
    numMuts = length(data$label)
  }
  maxMuts = max(table(data[, groupBy]))
  
  ## plot boundaries
  xmin = 0
  #xmax = max(max(data$pos), max(regions$end)) * 3
  xmax = max(bmExons$cds_end/3)
  xlength = xmax - xmin
  ymin = 0
  ymax = 30  ## upper margin of the plot (too large: much empty space at the top of the plot; too small: mutation labels do not fit into the plot region)
  y0 = 6
  y1 = 7.5
  y2 = 10
  axis_steps = 100

  ## design parameters
  labels.cex = cex
  lwd = 1.3
  pch = 21
  
  ## initialize plot window
  plot(c(xmin-xlength*0.2, xmax+xlength*0.2), c(ymin,ymax), type="n", xaxt="n", yaxt="n", xlab="", ylab="", frame.plot=FALSE, xaxs="i", yaxs="i", main=title)
  axis(side=1, at=c(1, seq(xmin+axis_steps, xmax, axis_steps)), pos=y0*0.7, cex.axis=cex)

  ## text size information is important to correctly adjust the lines and points to the plotted text
  textHeight = max(strheight(letters, cex=labels.cex))
  textWidth =  max(strwidth(letters, cex=labels.cex))
  textWidth = textWidth + (textWidth * 0.2)

  ## draw annotated regions
  polygon(x=c(xmin,xmin,xmax,xmax), y=c(y0,y1,y1,y0), col="black", border="black", density=0, lwd=lwd)
  text(x=xmin, y=(y0[1]+y1[1])/2, labels=bmGene$hgnc_symbol[1], pos=2, cex=cex)
  text(x=xmax, y=(y0[1]+y1[1])/2, labels=paste(bmExons$cds_length[1]/3,"aa",sep=""), pos=4, cex=cex)
  if(!any(is.na(regions))){
    for(i in 1:nrow(regions)){
      polygon(x=c(regions$start[i], regions$start[i], regions$end[i], regions$end[i]), y=c(y0,y1,y1,y0), border="black", col=as.character(regions$color)[i], lwd=lwd)
    }
  }
  ## Draw exon regions
  for(i in 1:nrow(bmExons)){
    start = bmExons$cds_start[i]/3
    end = bmExons$cds_end[i]/3
    polygon(x=c(start,start,end,end), y=c(y0*0.9,y0*0.7,y0*0.7,y0*0.9), col="lightgrey", lwd=lwd)
    if(textWidth < (end-start)){
      text((start+end)/2, y0*0.8, bmExons$rank[i], cex=cex)
    }
  }

  ## y-positions of labels and mutation marks
  if(horiz== FALSE){
    y3 = y2+((maxMuts+1)*textHeight)
  }else{
    y3 = y2+(2*textHeight)
  }
  
  ## set mutation labels (take care of overlapping labels using xLab function from TeachingDemos package)
  x = data$pos[mutGroup]
  labels = data$label[mutGroup]
  ## if the plotting device is too small, it may not be possible to align all mutation labels without overlaps (gives a warning in function spread.labs)
  ## thus, catch this warning and retry with smaller font size (90% of original)
  overlappingLabels = TRUE
  while(overlappingLabels == TRUE){
    overlappingLabels = FALSE
    catchUpdate = tryCatch(
             xLab <- spread.labs(x=x, mindiff=textWidth, maxiter=100000, min=xmin, max=xmax),
             warning = function(w){
               message("Plot region too small for all mutation labels: Reducing font size. It may be appropriate to change the width of the device and restart the plot.")
               overlappingLabels = TRUE
               labels.cex = labels.cex * 0.9
               textWidth = max(strwidth(letters, cex = labels.cex))
               textWidth = textWidth + (textWidth * 0.2)
               return(list(overlappingLabels, labels.cex, textWidth))
             }
             )
    ## overlaps
    if(is.list(catchUpdate)){
      overlappingLabels = catchUpdate[[1]]
      labels.cex = catchUpdate[[2]]
      textWidth = catchUpdate[[3]]
    }
    ## no overlaps
    if(is.vector(catchUpdate)){
      xLab = catchUpdate
    }
  }

  x = data.frame(pos_axis=x, pos_labels=xLab, color=data$color[mutGroup], group=data[mutGroup, groupBy], stringsAsFactors=FALSE)
  colnames(x) = c("pos_axis", "pos_labels", "color", groupBy)
  text(x$pos_labels, y3+(textHeight/2), labels, srt=90, cex=labels.cex, adj=c(0,0.5), offset=3)

  ## mark mutations (as colored points)
  pt.cex = labels.cex
  ## vertical alignment of mutation marks
  if(horiz == FALSE){
    for(i in 1:numMuts){
      lines(list(x=c(x$pos_axis[i],x$pos_axis[i],x$pos_labels[i],x$pos_labels[i]), y=c(y0,y1,y2,y3)), col=x$color[i], lwd=lwd)
      xm = x$pos_labels[i]
      ym = y2 + textHeight
      for(m in muts[[i]]){
        if(!is.na(m)){
          points(xm, ym, bg=as.character(subset(mutationInfo, mutation == m)$color), col="black", cex=pt.cex, pch=pch, )
          ym = ym + textHeight
        }
      }
  }
  ## horizontal alignment of mutation marks
  }else{
    ## draw single mutation lines
    x2 = x[!(x[, groupBy] %in% x[duplicated(x[groupBy]), groupBy]), ]
    for(i in 1:nrow(x2)){
      lines(list(x=c(x2$pos_axis[i],x2$pos_axis[i],x2$pos_labels[i],x2$pos_labels[i]), y=c(y0,y1,y2,y3)), col=x2$color[i], lwd=lwd)
    }
    ## draw lines for mutation groups (as "forks" like in dendrograms)
    groups = table(x[, groupBy])
    groups = groups[groups > 1]
    anyGroups = any(groups > 1)
    if(anyGroups == TRUE){
      groups = groups[groups > 1]
      x_group= x[x[, groupBy] %in% x[duplicated(x[groupBy]), groupBy], ]
      for(i in 1:length(groups)){
        x2 = x_group[x_group[, groupBy] == names(groups[i]), ]
        lines(x=c(min(x2$pos_labels), min(x2$pos_labels), max(x2$pos_labels), max(x2$pos_labels)), y=c(y3,y2,y2,y3), col=x2$color[1], lwd=lwd)
        #points(x=(max(x2$pos_labels)+min(x2$pos_labels))/2, y=y2, bg=x2$color[1], pch=25, cex=pt.cex*0.75)  ## "anchor point"
        polygon(x=c(min(x2$pos_labels), min(x2$pos_labels), max(x2$pos_labels), max(x2$pos_labels)), y=c(y2,y3,y3,y2), col=x2$color[1], border=NA, density=18, lwd=lwd*0.75)
        lines(list(x=c(x2$pos_axis[1],x2$pos_axis[1],(max(x2$pos_labels)+min(x2$pos_labels))/2), y=c(y0,y1,y2)), col=x2$color[1], lwd=lwd)
      }
    }
    ## draw mutations marks
    for(i in 1:nrow(x)){
      m = muts[[i]]
      xm = x$pos_labels[i]
      ym = y2 + textHeight
      points(xm, ym, bg=as.character(subset(mutationInfo, mutation == m)$color), col="black", cex=pt.cex, pch=pch)
    }
  }
  
  ## draw legend for mutations and annotated regions
  if(legend == TRUE){
    legend(x="topleft", legend=mutationInfo$legend, pt.bg=as.character(mutationInfo$color), col="black", pch=pch, bty="n", cex=cex, pt.cex=cex)
    if(!any(is.na(regions))){  
      legend(x="bottomleft", legend=regions$name, fill=as.character(regions$color), bty="n", cex=cex, horiz=TRUE)
    }
  }

  return(bmReturn)
}

.annotatedVariants_to_dataFrame <- function(varAnnot, transcript){

  ## convert input of class AnnotatedVariants
  data = annotatedVariants(varAnnot)
  mutNames = names(data)
  label = vector(mode="character")
  pos = vector(mode="integer")
  type = vector(mode="character")
  m_count = 1
#  aminos = getAminoAbbr()  ## maybe offer the possibility to label amino acid changes
  for(mn in mutNames){
    
    m = data[[mn]]$exons
    
    ## see if exon data available
    if(nrow(m) > 0){
    
      ## see if variant lies on the given transcript
      m = subset(m, m$ensembl_transcript_id == transcript)
      if(nrow(m) > 0){
        
        for(i in 1:nrow(m)){
          
          ## draw markers for variants that lie in a coding area
          if(m[i,"coding"] == TRUE){
            if(length(grep("-", m[i, "codonMut"])) == 0){
              pos[m_count] = m[i, "numCodonStart"]             
              if(m[i,"AminoRef"] != "*" & m[i,"AminoMut"] != "*"){                
                ## missense point mutation
                if(m[i,"AminoRef"] !=  m[i,"AminoMut"]){
                  type[m_count]= "M"
 #                 label[m_count] = paste(aminos[m[i,"AminoRef"]], aminos[m[i,"AminoMut"]], sep=" -> ")
                }else{
                  ## silent point mutation
                  type[m_count] = "S"
 #                 label[m_count] = aminos[m[i,"AminoRef"]]
                }                
                ## nonsense point mut
              } else {
                if(m[i,"AminoRef"] != "*" & m[i,"AminoMut"] == "*"){
                  type[m_count] = "N"
 #                 label[m_count] = paste(aminos[m[i,"AminoRef"]], aminos[m[i,"AminoMut"]], sep=" -> ")
                }else
                warning("variant ", mn, " is of unknown type")                
              }             
              ## deletion
            } else {
              pos[m_count] = (m[i, "numCodonStart"] +
                              m[i, "numCodonEnd"]) / 2
              type[m_count] = "D"
 #             label[m_count] = paste(aminos[m[i,"AminoRef"]], m[i,"AminoMut"], sep=" -> ")
            }
          }
          m_count = m_count + 1
        }
      }
    }
  }
  
  df = data.frame(label=rep(" ", length(pos)), pos=pos, mutation=type, color=rep("black", length(pos)), stringsAsFactors=FALSE)
  mutationInfo = data.frame(mutation=c("M","N","S","D"), legend=c("Missense","Nonsense","Silent","Deletion"), color=c("cadetblue4","darkred","darkgoldenrod1","blue"), stringsAsFactors=FALSE)
  return(list(df, mutationInfo))
}

.getEnsemblInfo <- function(gene, transcript){
  ## retrieve exon and gene information from Ensembl (if Ensembl gene-id is provided)
  ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl", host="http://grch37.ensembl.org")
  bmAttributes = c("ensembl_gene_id", "ensembl_transcript_id", "rank", "cds_start", "cds_end", "cds_length")
  bmExons = getBM(attributes=bmAttributes, filters="ensembl_gene_id", values=gene, mart=ensembl)
  bmExons = subset(bmExons, !apply(is.na(bmExons), 1, any))
  bmReturn = bmExons
  if(is.na(transcript)){
    ## take largest transcript if no ensembl transcript id is given
    transcript = subset(bmExons, cds_length == max(bmExons$cds_length)[1])$ensembl_transcript_id[1]
    warning("No Ensembl transcript-id given. Will choose largest transcript. It is recommended to specify the transcript via the \"transcript\" parameter to ensure a most appropriate illustration.")
  }
  bmExons = subset(bmExons, ensembl_transcript_id == transcript)
  bmAttributes = c("ensembl_gene_id", "hgnc_symbol")
  bmGene = getBM(attributes=bmAttributes, filters="ensembl_gene_id", values=gene, mart=ensembl)
  return(list(gene, transcript, bmExons, bmGene, bmReturn))
}




## Method for the R453Toolbox class "AnnotatedVariants"
setMethod("plotVariants",
          signature=c(data="AnnotatedVariants", "character"),
          function(data, gene, transcript=NA, regions, horiz=FALSE, cex=1, title="", legend=TRUE){
            
            #library(TeachingDemos)
            #library(biomaRt)

            ## check integrity of transcript
            if(!(is.character(transcript) | is.na(transcript))){
                stop("Invalid argument: transcript has to be an Ensembl transcript-id of type character")
            }
            ## check integrity of regions
            if(missing(regions)){
              regions = NA
            }else{
              if(!is.data.frame(regions))
                stop("Invalid argument: regions has to be of type data.frame")
              if(!all( c("name", "start", "end") %in% colnames(regions)))
                stop("Wrong columns in regions: please provide data as a data.frame with columns \"name\", \"start\", \"end\" and \"color\". See ?plotVariants for details.")
              if(!("color" %in% colnames(regions))){
                warning("No colors specified for given regions. Assigning colors automatically. If you find them inappropriate, please add a column \"color\" to your region data.frame. See parameter \"regions\" in ?plotVariants for details.")
                regions$color = rainbow(nrow(regions))
              }
            }
            ## check integrity of horiz, title, legend, cex
            if(!is.logical(horiz))
              stop("Invalid argument: horiz has to be TRUE or FALSE")
            if(!is.character(title))
              stop("Invalid argument: title has to be of type character")
            if(!is.logical(legend))
              stop("Invalid argument: legend has to be TRUE or FALSE")
            if(!is.numeric(cex))
              stop("Invalid argument: cex has to be a numeric value")

            ## get Ensembl info for given gene / transcript
            bmInfo = .getEnsemblInfo(gene, transcript)
            
            df = .annotatedVariants_to_dataFrame(data, bmInfo[[2]])
            data = df[[1]]
            mutationInfo = df[[2]]
            
            if(nrow(data) >0){
              .plotVariants(data=data, gene=bmInfo[[1]], transcript=bmInfo[[2]], bmExons=bmInfo[[3]], bmGene=bmInfo[[4]], bmReturn=bmInfo[[5]], regions=regions, mutationInfo=mutationInfo, groupBy="pos", cex=cex, horiz=horiz, title=title, legend=legend)
            }else{
              stop("No annotation for given transcript. Stopping here. I'm sorry.")
            }
            
          })


## Method for data in data.frame format
setMethod("plotVariants",
          signature=c(data="data.frame", gene="character"),
          function(data, gene, transcript=NA, regions, mutationInfo, groupBy="pos", horiz=FALSE, cex=1, title="", legend=TRUE){

            #library(TeachingDemos)
            #library(biomaRt)
        
            ## check integrity of data
            data = subset(data, !apply(is.na(data), 1, any))
            if(!all(c("label", "pos", "mutation") %in% colnames(data)))
              stop("Wrong data columns: please provide data as a data.frame with columns \"label\", \"pos\", \"mutation\" (optional: \"color\"). See ?plotVariants for details.")
            if(!("color" %in% colnames(data)))
              data$color = "black"
            
            ## check integrity of transcript
              if(!(is.character(transcript) | is.na(transcript))){
                stop("Invalid argument: transcript has to be an Ensembl transcript-id of type character")
            }
            ## check integrity of regions
            if(missing(regions)){
              regions = NA
            }else{
              if(!is.data.frame(regions))
                stop("Invalid argument: regions has to be of type data.frame")
              if(!all( c("name", "start", "end") %in% colnames(regions)))
                stop("Wrong columns in regions: please provide data as a data.frame with columns \"name\", \"start\", \"end\" and \"color\". See ?plotVariants for details.")
              if(!("color" %in% colnames(regions))){
                warning("No colors specified for given regions. Assigning colors automatically. If you find them inappropriate, please add a column \"color\" to your region data.frame. See parameter \"regions\" in ?plotVariants for details.")
                regions$color = rainbow(nrow(regions))
              }
            }
            ## check integrity of mutationInfo
            if(missing(mutationInfo)){
              warning("Legend and colors for all mutation types were generated automatically. It is suggested to specify your own mutation info. See parameter \"mutationInfo\" in ?plotVariants for details.")
              muts = levels(factor(data$mutation))
              col = rainbow(length(muts))
              mutationInfo = data.frame(mutation=muts, legend=muts, color=col, stringsAsFactors=FALSE)              
            }else{
              if(!is.data.frame(mutationInfo))
                stop("Invalid argument: mutationInfo has to be a data frame")
              if(!all(c("mutation", "legend") %in% colnames(mutationInfo)))
                stop("Wrong columns in mutationInfo: please provide data as a data.frame with columns \"mutation\", \"legend\" and \"color\". See ?plotVariants for details.")
              if(!("color" %in% colnames(mutationInfo))){
                warning("No colors specified for given mutation types. Assigning colors automatically. If you find them inappropriate, please add a column \"color\" to your mutationInfo data.frame. See parameter \"mutationInfo\" in ?plotVariants for details.")
                mutationInfo$color = rainbow(nrow(mutationInfo))
              }
            }
            ## check integrity of groupBy
            if(!is.character(groupBy)){
              stop("Invalid argument: groupBy has to be of type character")
            }else{
              if(!(groupBy %in% colnames(data)))
                stop("Column not found: When grouping mutations by a column other than \"pos\", make sure the column has the same name in your data")
            }
            ## check integrity of horiz, title, legend, cex
            if(!is.logical(horiz))
              stop("Invalid argument: horiz has to be TRUE or FALSE")
            if(!is.character(title))
              stop("Invalid argument: title has to be of type character")
            if(!is.logical(legend))
              stop("Invalid argument: legend has to be TRUE or FALSE")
            if(!is.numeric(cex))
              stop("Invalid argument: cex has to be a numeric value")

            ## get Ensembl info for given gene / transcript
            bmInfo = .getEnsemblInfo(gene, transcript)
            
            return(.plotVariants(data=data, gene=bmInfo[[1]], transcript=bmInfo[[2]], bmExons=bmInfo[[3]], bmGene=bmInfo[[4]], bmReturn=bmInfo[[5]], regions=regions, mutationInfo=mutationInfo, groupBy=groupBy, horiz=horiz, cex=cex, title=title, legend=legend))
          })
