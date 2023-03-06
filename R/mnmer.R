#' mnmer
#' @aliases cmnmer mnmer MNmer
#' @description Generates the feature matrix using conditional probability. As default, all sequences with N+IUPAC content higher 10% than are removed.
#'
#' @param FASTAfile FASTA file \code{FASTAfile}
#' @param k Int value of k \code{k}
#' @param m Int value of m \code{m}
#' @param pni Cutoff percentage of N+IUPAC bases
#'
#' @return Outputs a dataframe 
#'
#' @keywords (m,n)-mer k-mers featureMatrix
#'
#' @export
#' 
#' @examples
#' dir <-system.file("extdata", package="mnmer")
#' mosquito <- mnmer(file.path(dir, "mosquito_vir.fasta.gz"),2,0)
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

    fasfile <- .Call ("readrandFASTA", size, pni)


}

mnmer <- function (seqs, k, m) 
{
    k <- as.integer(k)
    m <- as.integer(m)
    if (typeof(k) !="integer" || typeof(m)!="integer" || k == 0)
        stop("ERROR: parameters k and m must be integer type.")
    
    if (k < m)
        stop ("k must be greater or equal to m")

    seqid <- c()
    tab <- data.frame()

    for (i in 1:length(seqs)){
        seqid <- c(seqid,names(seqs)[i])
        sqs <- as.character(seqs[[i]])
        ctab <- .Call("cmnmer", sqs, k, m)
        ts <- read.csv(text = ctab, header = TRUE)
        tab <- rbind (tab,ts)
    }

    tab <- cbind (seqid,tab)

    return(tab)
}


