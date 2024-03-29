% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/itempool-class-methods.R,
%   R/response-class-methods.R, R/response_set-class-methods.R,
%   R/testlet-class-methods.R
\name{length,Itempool-method}
\alias{length,Itempool-method}
\alias{length,Response-method}
\alias{length,Response_set-method}
\alias{length,Testlet-method}
\title{Find the length of an \code{\link{Itempool-class}} object}
\usage{
\S4method{length}{Itempool}(x)

\S4method{length}{Response}(x)

\S4method{length}{Response_set}(x)

\S4method{length}{Testlet}(x)
}
\arguments{
\item{x}{an \code{\link{Response_set-class}} object}
}
\description{
Find the length of an \code{\link{Itempool-class}} object

Find the length of an \code{\link{Response-class}} object

Find the length of a \code{\link{Response_set-class}} object

Find the length of a \code{\link{Testlet-class}} object
}
\examples{
r <- response(sample(0:1, 22, TRUE))
length(r)

}
\author{
Emre Gonulates
}
