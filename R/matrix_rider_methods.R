setMethod("getSeqOccupancy", signature=c(pwm="PFMatrix", cutoff="numeric"),
   function(pwm, cutoff) {
      c_res <- .Call("get_occupancy", pwm, cutoff)
      print(c_res)
      return(c_res)
   }
)
