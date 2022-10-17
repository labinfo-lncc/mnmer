#' mnmer
#' 
#' @description Generates the feature matrix using conditional probability. 
#'
#' @param inputParameter1 FASTA file \code{inputParameter1}
#' @param inputParameter2 Int value of k \code{inputParameter2}
#' @param inputParameter3 Int value of m \code{inputParameter3}
#'
#' @return Outputs a dataframe 
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
        stop("Fix the parameters")}
    ctab <- .Call("cmnmer", file, k, m)
    tab <- read.csv(text = ctab, header = TRUE)
    return(tab)
}


