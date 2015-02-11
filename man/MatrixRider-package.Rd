\name{MatrixRider-package}
\alias{MatrixRider-package}
\alias{MatrixRider}
\docType{package}
\title{
Calculate total affinity and occupancies for binding site matrices on a given sequence
}
\description{
Calculates a single number for a whole sequence that reflects the propensity of
a DNA binding protein to interact with it. The DNA binding protein has to be described
with a PFM matrix, for example gotten from Jaspar.
}
\details{
\tabular{ll}{
   Package: \tab MatrixRider\cr
   Type: \tab Package\cr
   Version: \tab 0.99.1 \cr
   Date: \tab 2015-23-01\cr
   License: \tab GPL-3\cr
   }
}

\author{
Elena Grassi <elena.grassi@unito.it>
}

\keyword{MatrixRider}
\examples{
      library(JASPAR2014)
      library(TFBSTools)
      library(Biostrings)
      pfm <- getMatrixByID(JASPAR2014,"MA0004.1")
      ## The following sequence has a single perfect match 
      ## thus it gives the same results with all cutoff values.
      sequence <- DNAString("CACGTG")
      getSeqOccupancy(sequence, pfm, 0.1)
      getSeqOccupancy(sequence, pfm, 1)
}
