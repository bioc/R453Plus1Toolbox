.calculateTiTv <- function(object) {

  if (!is.element("variantBase", names(fData(object))) |
      !is.element("referenceBases", names(fData(object)))) {
    stop("Column variantBase or referenceBases is missing.")
  }

  df = fData(object)
  indSNP = nchar(df$referenceBases) == 1 & nchar(df$variantBase) == 1
  df = df[indSNP, ]

  subMat = matrix(NA, nrow=4, ncol=4,
      dimnames=list(c("A", "G", "C", "T"), c("A", "G", "C", "T")))
  subMat["A", "G"] = sum(df$referenceBases == "A" & df$variantBase == "G")
  subMat["A", "C"] = sum(df$referenceBases == "A" & df$variantBase == "C")
  subMat["A", "T"] = sum(df$referenceBases == "A" & df$variantBase == "T")
  subMat["G", "A"] = sum(df$referenceBases == "G" & df$variantBase == "A")
  subMat["G", "C"] = sum(df$referenceBases == "G" & df$variantBase == "C")
  subMat["G", "T"] = sum(df$referenceBases == "G" & df$variantBase == "T")
  subMat["C", "A"] = sum(df$referenceBases == "C" & df$variantBase == "A")
  subMat["C", "G"] = sum(df$referenceBases == "C" & df$variantBase == "G")
  subMat["C", "T"] = sum(df$referenceBases == "C" & df$variantBase == "T")
  subMat["T", "A"] = sum(df$referenceBases == "T" & df$variantBase == "A")
  subMat["T", "G"] = sum(df$referenceBases == "T" & df$variantBase == "G")
  subMat["T", "C"] = sum(df$referenceBases == "T" & df$variantBase == "C")
    
  # A transition is a A <-> G or a C <-> T substitution
  # http://www.mun.ca/biology/scarr/Transitions_vs_Transversions.html
  # there are twice as many transversion
  # expected values: http://www.broadinstitute.org/gsa/wiki/index.php/VariantEval
  ti = sum(subMat["A", "G"] + subMat["G", "A"]
      + subMat["C", "T"] + subMat["T", "C"])
  tv = sum(subMat, na.rm=TRUE) - ti

  return(list(subMat=subMat, TiTv= ti/tv))
}


setMethod("calculateTiTv", signature=signature(object="AVASet"),
    .calculateTiTv)

setMethod("calculateTiTv", signature=signature(object="MapperSet"),
    .calculateTiTv)




