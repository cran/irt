% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resp_lik.R
\name{resp_lik}
\alias{resp_lik}
\alias{resp_lik,Item-method}
\alias{resp_lik,Itempool-method}
\alias{resp_lik,Testlet-method}
\title{Likelihood of a response string}
\usage{
resp_lik(ip, resp, theta)

\S4method{resp_lik}{Item}(ip, resp, theta)

\S4method{resp_lik}{Itempool}(ip, resp, theta)

\S4method{resp_lik}{Testlet}(ip, resp, theta)
}
\arguments{
\item{ip}{An \code{\link{Item-class}}, \code{\link{Itempool-class}} or a
\code{\link{Testlet-class}} object.}

\item{resp}{A vector of item responses.}

\item{theta}{An vector containing ability parameters.}
}
\value{
A matrix of likelihood(s)
}
\description{
\code{resp_lik} returns the likelihood of a response string
for given items and ability.
}
\examples{
item <- generate_item(model = "3PL")
theta <- rnorm(6)
resp <- sim_resp(ip = item, theta = theta, prop_missing = .1)
resp_lik(ip = item, resp = resp, theta = theta)

item <- generate_item(model = "GRM")
resp <- sim_resp(ip = item, theta = theta, prop_missing = .1)
resp_lik(ip = item, resp = resp, theta = theta)
ip <- generate_ip(model = "3PL")
theta <- rnorm(6)
resp <- sim_resp(ip = ip, theta = theta, prop_missing = .1)
resp_lik(ip = ip, resp = resp, theta = theta)

ip <- generate_ip(model = "GRM")
resp <- sim_resp(ip = ip, theta = theta, prop_missing = .1)
resp_lik(ip = ip, resp = resp, theta = theta)

}
\author{
Emre Gonulates
}
