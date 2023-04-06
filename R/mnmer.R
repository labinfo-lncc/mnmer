#' mnmer
#' @aliases cmnmer mnmer MNmer
#' @description Generates the feature matrix using conditional probability. As default, all sequences with N+IUPAC content higher 10% than are removed.
#'
#' @param seqs DNAStringSet object \code{seqs}
#' @param m Int value of m \code{m}
#' @param n Int value of n \code{n}
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

mnmer <- function (seqs, m, n) 
{
    m <- as.integer(m)
    n <- as.integer(n)

    if (typeof(m) !="integer" || typeof(n)!="integer" || m == 0)
        stop("ERROR: parameters m and n must be integer types and m must be different zero.")
    
    seqid <- c()
    tab <- data.frame()
    pb = txtProgressBar(min = 0, max = length(seqs), initial = 0) 

    for (i in 1:length(seqs)){
        seqid <- c(seqid,names(seqs)[i])
        sqs <- as.character(seqs[[i]])
        ctab <- .Call("cmnmer", sqs, m, n)
        ts <- read.csv(text = ctab, header = TRUE)
        tab <- rbind (tab,ts)
        setTxtProgressBar(pb,i)
    }

    tab <- cbind (seqid,tab)

    return(tab)
}



readNumFASTA <- function (FASTAfile, size=0, rand=FALSE, pni=0.02)
{
    if (!file.exists(FASTAfile)) 
        stop(paste0("File: ", FASTAfile, " does not found!"))
    if (!rand) {
        seqs <- readDNAStringSet(FASTAfile)
        seqid <- c()
        tab <- data.frame()
        maxsize <- ifelse(size == 0, length(seqs), size)
        dsqs <- DNAStringSet()
        for (i in 1:maxsize) {
            ni <- alphabetFrequency(seqs[i], baseOnly = TRUE, 
                as.prob = TRUE)
            if (ni[5] < pni) {
                dsqs <- c(dsqs, seqs[i])
            }
            else {
                message("Warning: ", names(seqs)[i], " has a proportion of N + IUPAC bases = ", 
                  ni[5])
            }
        }
        return(dsqs)
    }
    fasstr <- .Call("readrandFASTA", FASTAfile, size, pni)
    fastab <- read.table(text = fasstr, header = FALSE, sep = "\t")
    dsqs <- DNAStringSet(fastab[, 2])
    names(dsqs) <- fastab[, 1]
    return(dsqs)
}




