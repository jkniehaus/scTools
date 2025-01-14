% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seuratHelpR.R
\name{StackedVlnPlot}
\alias{StackedVlnPlot}
\title{Stacked Violin}
\usage{
StackedVlnPlot(
  obj,
  features,
  pt.size = 0,
  plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
  ...
)
}
\arguments{
\item{obj}{Seurat object}

\item{features}{gene names}

\item{plot.margin}{Plot Margins}
}
\description{
This function creates a stacked violin plot from the seurat VlnPlot function
}
\examples{
StackedVlnPlot()
}
\keyword{stacked}
\keyword{violin}
