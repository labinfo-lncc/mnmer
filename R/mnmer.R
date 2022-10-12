#' The mn function generates the feature matrix using conditional probability.
#' 
#' This function can generate both k-mers and (m,n)-mers. It calls the cmnmer script. 
#'
#' @param inputParameter1 FASTA file \code{inputParameter1}
#' @param inputParameter2 Value of k for k-mer generation. Needs to be different from zero. \code{inputParameter2}
#' @param inputParameter3 Value of m for (m,n)-mer generation in the format of (m, k-m). In case of k-mer generation, m should be zero as (0,k). \code{inputParameter3}
#'
#' @return output A description of the object the function outputs 
#'
#' @keywords (m,n)-mer k-mers featureMatrix
#'
#' @export
#' 
#' @examples
#' features <- mnmer(sequences.fasta,2,0)
#' features <- mnmer(sequences.fasta,2,1)
#' 
mnmer <- function (file, k, m) {
  if (!file.exists(file)) {
    stop("The file does not exists")}
  
  k <- as.integer(k)
  m <- as.integer(m)
  
  if (typeof(k) !="integer" || typeof(m)!="integer" || k =="0" || m > k){
    stop("Fix the parameters")
  }
  
    ctab <- .Call("cmnmer", file, k, m)
    tab <- read.csv(text = ctab, header = T)
    return(tab)
}


