\name{getSeqOccupancy}
   \alias{getSeqOccupancy}
   \alias{getSeqOccupancy-methods}
   \alias{getSeqOccupancy,DNAString,PFMatrix,numeric-method}
   \alias{getSeqOccupancy,DNAString,PFMatrixList,numeric-method}
   \title{TODO}
   \description{
       TODO
      }
   \usage{
      getSeqOccupancy(sequence, pfm, cutoff)
   }
   \arguments{
      \item{sequence}{
          A \linkS4class{DNAString} object TODO
         }
      \item{pfm}{
         A \linkS4class{PFMatrix} or a \linkS4class{PFMatrixList} object TODO
         }
      \item{cutoff}{
         TODO
         }
      }
   \value{
      TODO
   }
   \examples{
      library(JASPAR2014)
      library(TFBSTools)
      library(Biostrings)
      pfm <- getMatrixByID(JASPAR2014,"MA0004.1")
      sequence <- DNAString("GAATTCTCTCTTGCACGTGTAGTCTCTTGACAAAATG")
      getSeqOccupancy(sequence, pfm, 0.1)
   }
\keyword{getSeqOccupancy}
