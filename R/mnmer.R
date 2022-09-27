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


