setMethod("getSeqOccupancy", signature=c(sequence="DNAString", pfm="PFMatrix", cutoff="numeric"),
   function(sequence, pfm, cutoff) {
   cRes <- .Call("get_occupancy", pfm, cutoff, sequence)
   return(cRes)
   }
)

setMethod("getSeqOccupancy", signature=c(sequence="DNAString", pfm="PFMatrixList", cutoff="numeric"),
   function(sequence, pfm, cutoff){
      ansList = vapply(pfm, FUN=getSeqOccupancy, FUN.VALUE=1, cutoff=cutoff, sequence=sequence)
      return(ansList)
   }
)