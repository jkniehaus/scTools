% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seuratHelpR.R
\name{upset_plot}
\alias{upset_plot}
\title{upset plot}
\usage{
upset_plot(object, genes, metacol = NULL, metavar = NULL)
}
\arguments{
\item{object}{Seurat object}

\item{genes}{Gene names}

\item{metacol}{Column in metadata to subset cells by}

\item{metavar}{group in metadata column to subset cells by}
}
\description{
This function generates upset plots with UpSetR package
}
\examples{
upset_plot()
}
\keyword{upset}
