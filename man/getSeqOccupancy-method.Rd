\name{getSeqOccupancy}
\alias{getSeqOccupancy}
\alias{getSeqOccupancy-methods}
\alias{getSeqOccupancy,DNAString,PFMatrix,numeric-method}
\alias{getSeqOccupancy,DNAString,PFMatrixList,numeric-method}
\title{Computes the total affinity or the occupancy at a given cutoff}
\description{
    The affinity/occupancy is calculated on the given DNAString for the 
    given PFMatrix (or all the ones in the PFMatrix list given).
}
\usage{
    getSeqOccupancy(sequence, pfm, cutoff)
}
\arguments{
    \item{sequence}{
        A \linkS4class{DNAString} object with the sequence for which affinity
        will be computed.
    }
\item{pfm}{
    A \linkS4class{PFMatrix} or a \linkS4class{PFMatrixList} object 
    with the matrixes whose affinity will be calculated.
    The background (\link{bg,XMatrix-method}) of the given pfm is used 
    to perform affinity calculations.
    }
    \item{cutoff}{
    numeric(1); between 0 and 1 (included): 0 corresponds to total affinity 
    (i.e. summing all the affinities) while 1 to
    summing only values corresponding to the perfect match 
    for a given \linkS4class{PFMatrix}.
    See vignette for details on how scores are calculated. If
    \pkg{MatrixRider} is installed, open the vignette with
    \code{vignette("MatrixRider")}.
    }
}
\value{
numeric; the resulting total affinity calculated on the given fasta.
If a \linkS4class{PFMatrixList} has been passed then a named 
numeric vector with the affinities for all the PFMs.
The vignette has all the details on the calculations 
(such as PFM to PWM conversion and pseudocounts).
}
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

    pfm2 <- getMatrixByID(JASPAR2014,"MA0005.1")
    pfms <- PFMatrixList(pfm, pfm2)
    names(pfms) <- c(name(pfm), name(pfm2))
## This calculates total affinity for both the PFMatrixes.
    getSeqOccupancy(sequence, pfms, 0)
}
\keyword{getSeqOccupancy}
