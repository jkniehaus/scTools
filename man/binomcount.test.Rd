% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seuratHelpR.R
\name{binomcount.test}
\alias{binomcount.test}
\title{Binomial test}
\usage{
binomcount.test(object, cells.1, cells.2, effect.size)
}
\arguments{
\item{object}{Seurat object}

\item{cells.1}{Group 1; a working group within seurat object levels/Idents}

\item{cells.2}{Group 2; if NULL = all other cells}

\item{effect.size}{Filter DEGs by effect size}
}
\description{
This function tests the probability of a difference in the proportion of non-zero values b/w groups; used by other functions. Adopted from Shekhar 2016
}
\examples{
binomcount.test()
}
\keyword{binom}
