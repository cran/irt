% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/item_fit.R
\name{item_fit}
\alias{item_fit}
\title{Calculate item-fit indices}
\usage{
item_fit(ip, resp, theta, type = "Q3")
}
\arguments{
\item{ip}{An \code{\link{Itempool-class}} object.}

\item{resp}{A vector of item responses.}

\item{theta}{An vector containing ability parameters.}

\item{type}{The type of the item-fit index.}
}
\value{
A vector of item-fit index values.
}
\description{
\code{item_fit} calculates the fit of an item to a given psychometric model.
}
\author{
Emre Gonulates
}
