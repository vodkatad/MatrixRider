setMethod("getSeqOccupancy", signature=c(sequence="DNAString", pfm="PFMatrix", cutoff="numeric"),
   function(sequence, pfm, cutoff) {
      if (cutoff < 0 || cutoff > 1) {
         stop("Wrong argument to getSeqOccupancy, 'cutoff' has to be between 0 and 1! (0 <= cutoff <= 1)")
      }
      # Probably testing this here is better for the user but bad in terms of efficiency. XXX FIXME DO at the C level.
      if (! (all(vapply(seq(1, length(sequence)), FUN=function(x){letter(sequence,x) %in% c("A","C","G","T", "N")}, FUN.VALUE=FALSE)))) {
         stop("Wrong argument to getSeqOccupancy, 'sequence' must be based on a restricted alphabet with only 'A','C','G','T' and 'N'")
      }
      cRes <- .Call("get_occupancy", pfm, cutoff, sequence)
      return(cRes)
   }
)

setMethod("getSeqOccupancy", signature=c(sequence="DNAString", pfm="PFMatrixList", cutoff="numeric"),
   function(sequence, pfm, cutoff){
      if (cutoff < 0 || cutoff > 1) {
         stop("Wrong argument to getSeqOccupancy, 'cutoff' has to be between 0 and 1! (0 <= cutoff <= 1)")
      }
      if (! (all(vapply(seq(1, length(sequence)), FUN=function(x){letter(sequence,x) %in% c("A","C","G","T", "N")}, FUN.VALUE=FALSE)))) {
         stop("Wrong argument to getSeqOccupancy, 'sequence' must be based on a restricted alphabet with only 'A','C','G','T' and 'N'")
      }
      ansList = vapply(pfm, FUN=getSeqOccupancy, FUN.VALUE=1, cutoff=cutoff, sequence=sequence)
      return(ansList)
   }
)