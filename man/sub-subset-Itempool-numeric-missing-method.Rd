% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Itempool-class-methods.R
\name{[[<-,Itempool,numeric,missing-method}
\alias{[[<-,Itempool,numeric,missing-method}
\title{Set the elements of an \code{Itempool} objects.}
\usage{
\S4method{[[}{Itempool,numeric,missing}(x, i, j) <- value
}
\arguments{
\item{x}{\code{\link{Item-class}} object from which to extract element(s) or
in which to replace element(s).}

\item{i}{indices specifying elements to extract or replace.}

\item{j}{This will not be used in \code{\link{Itempool-class}} objects.}

\item{value}{An \code{\link{Item-class}} object.}

\item{...}{Additional parameters to be passed to the function.}
}
\value{
An \code{\link{Itempool-class}} object with elements from
  \code{\link{Item-class}} class.
}
\description{
Set the elements of an \code{Itempool} objects.
}
\examples{
item1 <- item(a = 1.12, b = -2.1, c = 0.28)
item2 <- item(a = 2, b = 3.2, c = 0.21)
ip <- c(item1, item2)
item3 <- item(a = 1, b = -.2, c = 0.4)
ip[[2]] <- item3
}
\author{
Emre Gonulates
}
\keyword{internal}
