#' mnmer
#' @aliases cmnmer mnmer MNmer
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
#' dir <-system.file("extdata", package="mnmer")
#' mosquito <- mnmer(file.path(dir, "mosquito_vir.fasta"),2,0)
#' 
mnmer <- function (inputParameter1, inputParameter2, inputParameter3) {
    if (!file.exists(inputParameter1)) {
        stop("The file does not exists")}
    inputParameter2 <- as.integer(inputParameter2)
    inputParameter3 <- as.integer(inputParameter3)
    if (typeof(inputParameter2) !="integer" || typeof(inputParameter3)!="integer" || inputParameter2 =="0" || inputParameter3 > inputParameter2){
        stop("Fix the parameters")}
    ctab <- .Call("cmnmer", inputParameter1, inputParameter2, inputParameter3)
    tab <- read.csv(text = ctab, header = TRUE)
    return(tab)
}


