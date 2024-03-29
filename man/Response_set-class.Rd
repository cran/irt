% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/response_set-class.R
\docType{class}
\name{Response_set-class}
\alias{Response_set-class}
\title{An S4 class representing responses of a set of examinees}
\description{
An S4 class representing responses of a set of examinees
}
\section{Slots}{

\describe{
\item{\code{response_list}}{A list of \code{\link{Response-class}} objects. If
the \code{examinee_id} slots of \code{\link{Response-class}} objects are
not \code{NULL}, there cannot be duplicates.}

\item{\code{item_id}}{A character vector of \code{Item} ID's in the
\code{\link{Response-class}} objects. The order of this item_id will be
used when converting \code{\link{Response_set-class}} objects to a matrix.}

\item{\code{testlet_id}}{A character vector of \code{Testlet} ID's in the
\code{\link{Response-class}} objects.}

\item{\code{misc}}{This slot will hold any other information about the response
set.}
}}

\author{
Emre Gonulates
}
