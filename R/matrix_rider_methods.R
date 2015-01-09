setMethod("getSeqOccupancy", signature=c(pwm="PFMatrix", cutoff="numeric"),
   function(pwm, cutoff, subject) {
      # DNAString, XStringViews or MaskedDNAString
      c_res <- .Call("get_occupancy", pwm, cutoff, subject)
      print(c_res)
      return(c_res)
   }
)
