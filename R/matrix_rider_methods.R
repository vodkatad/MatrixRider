setMethod("getSeqOccupancy", signature=c(sequence="DNAString", 
                                        pfm="PFMatrix",
                                        cutoff="numeric"),
    function(sequence, pfm, cutoff) {
        if (cutoff < 0 || cutoff > 1) {
            stop(paste0("Wrong argument to getSeqOccupancy, ",
                "'cutoff' has to be between 0 and 1! (0 <= cutoff <= 1)"))
        }
        # Probably testing this here is better for the user
        # but bad in terms of efficiency.
        # XXX FIXME DO at the C level.
        cRes <- .Call("get_occupancy", pfm, cutoff, sequence)
        return(cRes)
    }
)

setMethod("getSeqOccupancy", signature=c(sequence="DNAString",
                                        pfm="PFMatrixList",
                                        cutoff="numeric"),
    function(sequence, pfm, cutoff) {
        if (cutoff < 0 || cutoff > 1) {
            stop(paste0("Wrong argument to getSeqOccupancy, ",
                "'cutoff' has to be between 0 and 1! (0 <= cutoff <= 1)"))
        }
        ansList = vapply(pfm, FUN=getSeqOccupancy, FUN.VALUE=1,
                        cutoff=cutoff, sequence=sequence)
        return(ansList)
    }
)