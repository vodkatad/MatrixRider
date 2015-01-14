#!/usr/bin/env Rscript
# Script to perform matrix_rider -w -1 -y analysis

checkReadable <- function(filename) {
   res <- file.access(names=filename, mode=4) == 0
   if (!res) {
      warning(paste(filename, "is not readable", sep=" "))
   }
   res
}

arguments <- matrix(c(
   'help', 'h', 0, "logical",
   'debug', 'd', 1, "character",
   'fasta' , 'f', 1, "character",
   'cutoff'  , 'c', 1, "numeric"
), ncol=4, byrow=T)

library(getopt)
opt <- getopt(arguments)

if (!is.null(opt$help)) {
   stop(getopt(arguments, command=get_Rscript_filename(), usage=TRUE))
}

if (is.null(opt$fasta)) {
   stop("Missing fasta [-f filename] parameter\n")
}

if (is.null(opt$cutoff)) {
   stop("Missing cutoff [-c X where 0 <= x <= 1] parameter\n")
}

library(matrix_rider)
library(JASPAR2014)
library(TFBSTools)

if (!checkReadable(opt$fasta)) {
   stop("The given fasta file does not exist or is not readable")  
}

if (opt$cutoff > 1 || opt$cutoff < 0) {
   stop("Wrong cutoff [-c X where 0 <= x <= 1] parameter\n")
}
   
opts<-list()
opts['collection']="CORE"
opts['matrixtype']="PFM"
mat <- getMatrixSet(JASPAR2014, opts)
fa <- readDNAStringSet(opt$fa)
res<-lapply(fa, getSeqOccupancy, pfm=mat, cutoff=0)
write.table(t(as.dataframe(res)), sep="\t", quote=FALSE)
         

if (!is.null(opt$debug)) {
   save.image(file=opt$debug)
}