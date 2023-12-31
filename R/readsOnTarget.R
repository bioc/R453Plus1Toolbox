.readsOnTarget <- function(alnReads, targetRegion) {

    result = list()
    for (i in 1:length(alnReads)) {
        rd = GRanges(IRanges(start=alnReads[[i]]$pos,
                             end=alnReads[[i]]$pos + alnReads[[i]]$qwidth - 1),
                             seqnames=alnReads[[i]]$rname)
      
        result[[i]] = overlapsAny(rd, targetRegion)
    }

    return(result)
}

setMethod("readsOnTarget",
    signature=signature(alnReads="list", targetRegion="GRanges"),
    .readsOnTarget)



# this version sets all duplicates TRUE as soon as one matches the target
#onTargetComplex <- function(samDF, chipDesignBed, trackPos=1) {
#    require("rtracklayer")
#
#    refBed=import.ucsc(chipDesignBed, subformat="bed")
#    irCD=refBed[[trackPos]]$ranges
#    spacesCD=gsub("chr", "", refBed[[trackPos]]$space)
#    rlCD=IRangesList(irCD)
#    names(rlCD)=spacesCD
#
#    du1=duplicated(samDF[,"name"], fromLast=FALSE)
#    du2=duplicated(samDF[,"name"], fromLast=TRUE)
#    du=du1 | du2
#    samu=samDF[!du,]
#    samd=samDF[du,]
#
#    iru=IRanges(start=samu[,"start"], end=samu[,"end"])
#    rdu=RangedData(iru, space=samu[,"chr"])
#    rlu=ranges(rdu)
#    mau=match(rlu, rlCD)
#    mau=!is.na(mau)
#    samu[,ncol(samu)+1]=mau
#    colnames(samu)[ncol(samu)]="onTarget"
#        
#    ird=IRanges(start=samd[,"start"], end=samd[,"end"])
#    rdd=RangedData(ird, space=samd[,"chr"])
#    rld=ranges(rdd)
#    mad=match(rld, rlCD)
#    mad=!is.na(mad)
#    samd[,ncol(samd)+1]=mad
#    colnames(samd)[ncol(samd)]="onTarget"
#
#    for(name in unique(samd[,"name"])){
#      actrows=rownames(samd[samd[,"name"] == name,])
#      samd[actrows,"onTarget"]=any(samd[actrows,ncol(samd)])
#    }
#
#    sam=rbind(samu, samd)
#    reads=nrow(sam)
#    on=sum(sam[,"onTarget"])
#    off=sum(!sam[,"onTarget"])
#    onperc=on/reads*100
#    offperc=off/reads*100
#}
