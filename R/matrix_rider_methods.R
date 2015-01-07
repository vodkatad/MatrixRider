setMethod("getSeqOccupancy", signature=(pwm="PWMatrix"),
   function(pwm, subject, seqname="Unknown", min.score=0.8) {
      c_res <- .Call("get_occupancy", pwm)
      print(c_res)
      return(c_res)
   }
)
