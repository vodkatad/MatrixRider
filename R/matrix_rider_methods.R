setMethod("getSeqOccupancy", signature=c(pwm="PFMatrix", cutoff="numeric", sequence="DNAString"),
   function(pwm, cutoff, sequence) {
      # DNAString, XStringViews or MaskedDNAString
      c_res <- .Call("get_occupancy", pwm, cutoff, sequence)
      print(c_res)
      return(c_res)
   }
)
