setMethod("getSeqOccupancy", signature=c(sequence="DNAString", pfm="PFMatrix", cutoff="numeric"),
   function(sequence, pfm, cutoff) {
   cRes <- .Call("get_occupancy", pfm, cutoff, sequence)
   return(cRes)
   }
)

setMethod("getSeqOccupancy", signature=c(sequence="DNAString", pfm="PFMatrixList", cutoff="numeric"),
   function(sequence, pfm, cutoff){
      ansList = lapply(pfm, getSeqOccupancy, cutoff=cutoff, sequence=sequence)
      # try vapply(cutoffs, FUN=function(x) {getSeqOccupancy(sequence=sequence, pfm=pfm, cutoff=x)}, FUN.VALUE=1)
      return(unlist(ansList))
   }
)