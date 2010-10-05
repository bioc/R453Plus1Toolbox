getAminoAbbr <- function() {
  aminos = c("Alanine", "Arginine", "Asparagine", "Aspartic acid", "Cysteine", "Glutamic acid", 
             "Glutamine", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine", 
             "Methionine", "Phenylalanine", "Proline", "Serine", "Threonine", "Tryptophan", 
             "Tyrosine", "Valine", "Stop")
  names(aminos) = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
  return(aminos)
}
