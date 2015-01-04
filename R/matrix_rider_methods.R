setMethod("getSeqOccupancy", signature=(pwm="integer"),
   function(pwm, subject, seqname="Unknown", min.score=0.8) {
      c_res <- .Call("get_occupancy", pwm)
      cat(c_res)
      cat("\n")
      return(c_res)
   }
)
