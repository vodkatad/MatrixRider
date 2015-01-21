setMethod("getSeqOccupancy", signature=c(sequence="DNAString", pfm="PFMatrix", cutoff="numeric"),
   function(sequence, pfm, cutoff) {
   # DNAString, XStringViews or MaskedDNAString
   c_res <- .Call("get_occupancy", pfm, cutoff, sequence)
   return(c_res)
   }
)

setMethod("getSeqOccupancy", signature=c(sequence="DNAString", pfm="PFMatrixList", cutoff="numeric"),
   function(sequence, pfm, cutoff){
      ans_list = lapply(pfm, getSeqOccupancy, cutoff=cutoff, sequence=sequence)
      #ans = do.call(SiteSetList, ans_list)
      return(unlist(ans_list))
   }
)