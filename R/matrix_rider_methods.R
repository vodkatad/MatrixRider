setMethod("getSeqOccupancy", signature=c(sequence="DNAString", pfm="PFMatrix", cutoff="numeric"),
   function(sequence, pfm, cutoff) {
   # DNAString, XStringViews or MaskedDNAString
   cRes <- .Call("get_occupancy", pfm, cutoff, sequence)
   return(cRes)
   }
)

setMethod("getSeqOccupancy", signature=c(sequence="DNAString", pfm="PFMatrixList", cutoff="numeric"),
   function(sequence, pfm, cutoff){
      ansList = lapply(pfm, getSeqOccupancy, cutoff=cutoff, sequence=sequence)
      return(unlist(ansList))
   }
)