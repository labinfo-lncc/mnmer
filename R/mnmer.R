mnmer <- function (file, k, m)
{
	ctab <- .Call ("cmnmer", file, k, m)
	tab <- read.csv (text=ctab, header=T)

	return (tab)
}
