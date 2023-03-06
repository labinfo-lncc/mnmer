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
mnmer <- function (FASTAfile, k, m, pni=0.1) {
    if (!file.exists(FASTAfile)) {
        stop("The file does not exists")
    }
    k <- as.integer(k)
    m <- as.integer(m)
    if (typeof(k) !="integer" || typeof(m)!="integer" || k =="0" || m > k){
        stop("Fix the parameters")
    }

    seqs <- readDNAStringSet (FASTAfile)
    seqid <- c()
    tab <- data.frame()

    for (i in 1:length(seqs)){
        ni <-  alphabetFrequency (seqs[i], baseOnly=TRUE, as.prob=TRUE)
        if (ni[5] < pni){ 
            seqid <- c(seqid,names(seqs)[i])
            sqs <- as.character(seqs[[i]])
            ctab <- .Call("cmnmer", sqs, k, m)
            ts <- read.csv(text = ctab, header = TRUE)
            tab <- rbind (tab,ts)
        }
        else {
           message("Warning: ", names(seqs)[i], " has a proportion of N + IUPAC bases = ", ni[5])
        }
    }

    tab <- cbind (seqid,tab)

    return(tab)
}


