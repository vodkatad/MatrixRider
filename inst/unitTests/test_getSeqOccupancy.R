# library(RUnit)
# source('inst/unitTests/test_getSeqOccupancy.R')
# test_getSeqOccupancy()
test_getSeqOccupancy_singlePWM_manyCutoffs <- function() {
   library(JASPAR2014)
   library(TFBSTools)
   library(Biostrings)
   pfm <- getMatrixByID(JASPAR2014,"MA0004.1")
   sequence <- DNAString("CACGTG")
   cutoffs <- seq(0,1,0.1)
   given <- vapply(cutoffs, FUN=function(x) {getSeqOccupancy(sequence=sequence, pfm=pfm, cutoff=x)}, FUN.VALUE=1)
   wanted <- rep(1470.946,11)
   checkEqualsNumeric(given, wanted, tolerance=1e-5)
}
