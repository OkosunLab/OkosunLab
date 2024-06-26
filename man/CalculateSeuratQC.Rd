% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleCellHelpers.R
\name{CalculateSeuratQC}
\alias{CalculateSeuratQC}
\title{CalculateSeuratQC}
\usage{
CalculateSeuratQC(seuratObj)
}
\arguments{
\item{seuratList}{a list of seurat objects}
}
\value{
a list of seurat objects
}
\description{
A function to Calculate QC metrics from a list  of seurat objects
Useful if you're stuck thinking of a palette
}
\examples{
CalculateSeuratQC()

CalculateSeuratQC(seuList)

}
\keyword{scRNAseq}
