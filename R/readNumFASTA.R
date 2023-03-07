#
 readNumFASTA
#' @aliases readNumFASTA readrandFASTA MNmer
#' @description Generates the feature matrix using conditional probability. As default, all sequences with N+IUPAC content higher 10% than are removed.
#'
#' @param FASTAfile FASTA file \code{FASTAfile}
#' @param size Int value of size \code{size}
#' @param rand Bool value of rand \code{rand}
#' @param pni Cutoff percentage of N+IUPAC bases
#'
#' @return Outputs a Biostring DNAStringSet
#'
#' @keywords FASTA Filter Random
#'
#' @export
#' 
#' @examples
#' dir <-system.file("extdata", package="mnmer")
#' mosquito <- readNumFASTA (file.path(dir, "mosquito_vir.fasta.gz"),size=0,rand=FALSE)
#' 
#' 

readNumFASTA <- function (FASTAfile, size=0, rand=FALSE, pni=0.02)
{
    if (!file.exists(FASTAfile))
        stop(paste0("File: ", FASTAfile, " does not found!"))
    
    if (!rand){
        seqs <- readDNAStringSet (FASTAfile)
        seqid <- c()
        tab <- data.frame()

        maxsize <- ifelse (size == 0, length(seqs), size)
        dsqs <- DNAStringSet()
        
        for (i in 1:maxsize){
            ni <-  alphabetFrequency (seqs[i], baseOnly=TRUE, as.prob=TRUE)

            if (ni[5] < pni){
                dsqs <- c(dsqs, seqs[i])
            }
            else {
                message("Warning: ", names(seqs)[i], " has a proportion of N + IUPAC bases = ", ni[5])
            }
        }    
        return(dsqs)
    }

    fasstr <- .Call ("readrandFASTA",FASTAfile,size, pni)
    fastab <- read.table(text = fasstr, header = FALSE, sep="\t")

    dsqs <- DNAStringSet (fastab[,2])
    names(dsqs) <- fastab[,1]

    return(dsqs)
}

