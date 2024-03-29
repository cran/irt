% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob.R
\name{prob}
\alias{prob}
\alias{prob,Item-method}
\alias{prob,Rasch-method}
\alias{prob,1PL-method}
\alias{prob,2PL-method}
\alias{prob,3PL-method}
\alias{prob,4PL-method}
\alias{prob,GRM-method}
\alias{prob,PCM-method}
\alias{prob,GPCM-method}
\alias{prob,GPCM2-method}
\alias{prob,Itempool-method}
\alias{prob,Testlet-method}
\alias{prob,numMatDfListChar-method}
\title{Calculate the probability of a correct response}
\usage{
prob(ip, theta, derivative = 0)

\S4method{prob}{Item}(ip, theta, derivative = 0)

\S4method{prob}{Rasch}(ip, theta, derivative = 0)

\S4method{prob}{1PL}(ip, theta, derivative = 0)

\S4method{prob}{2PL}(ip, theta, derivative = 0)

\S4method{prob}{3PL}(ip, theta, derivative = 0)

\S4method{prob}{4PL}(ip, theta, derivative = 0)

\S4method{prob}{GRM}(ip, theta, derivative = 0)

\S4method{prob}{PCM}(ip, theta, derivative = 0)

\S4method{prob}{GPCM}(ip, theta, derivative = 0)

\S4method{prob}{GPCM2}(ip, theta, derivative = 0)

\S4method{prob}{Itempool}(ip, theta, derivative = 0)

\S4method{prob}{Testlet}(ip, theta, derivative = 0)

\S4method{prob}{numMatDfListChar}(ip, theta, derivative = 0)
}
\arguments{
\item{ip}{An \code{\link{Item-class}}, or an \code{\link{Itempool-class}} or
\code{\link{Testlet-class}}  object containing the item parameters.}

\item{theta}{An object containing the ability parameters.}

\item{derivative}{Whether to calculate the first or second derivative of
probability of a response.
\describe{
  \item{\code{0}}{No derivative will be calculated. This is the default
    value.}
  \item{\code{1}}{Calculate the first derivative.}
  \item{\code{2}}{Calculate the second derivative.}
}}
}
\value{
Item probabilities at given theta will be returned.
}
\description{
\code{prob} Returns the probability of correct respond to an
  item or multiple items with given parameters for a given ability or
  abilities, i.e. \eqn{\theta}. For polytomous models, where there are
  multiple possible responses, probability of each response category will be
  returned.
}
\examples{
theta <- rnorm(1)
item1 <- generate_item(model = "Rasch")

# Probability of correct response
prob(item1, theta)

# First derivative of probability of correct response:
prob(item1, theta, derivative = 1)

# Second derivative of probability of correct response:
prob(item1, theta, derivative = 2)


# Multiple theta values
theta_n <- rnorm(5)

prob(item1, theta_n)
prob(item1, theta_n, derivative = 1)
prob(item1, theta_n, derivative = 2)

theta <- rnorm(1)
item1 <- generate_item(model = "1PL")

# Probability of correct response
prob(item1, theta)

# First derivative of probability of correct response:
prob(item1, theta, derivative = 1)

# Second derivative of probability of correct response:
prob(item1, theta, derivative = 2)


# Multiple theta values
theta_n <- rnorm(5)

prob(item1, theta_n)
prob(item1, theta_n, derivative = 1)
prob(item1, theta_n, derivative = 2)

theta <- rnorm(1)
item1 <- generate_item(model = "2PL")

# Probability of correct response
prob(item1, theta)

# First derivative of probability of correct response:
prob(item1, theta, derivative = 1)

# Second derivative of probability of correct response:
prob(item1, theta, derivative = 2)


# Multiple theta values
theta_n <- rnorm(5)

prob(item1, theta_n)
prob(item1, theta_n, derivative = 1)
prob(item1, theta_n, derivative = 2)

theta <- rnorm(1)
item1 <- generate_item(model = "3PL")

# Probability of correct response
prob(item1, theta)

# First derivative of probability of correct response:
prob(item1, theta, derivative = 1)

# Second derivative of probability of correct response:
prob(item1, theta, derivative = 2)


# Multiple theta values
theta_n <- rnorm(5)

prob(item1, theta_n)
prob(item1, theta_n, derivative = 1)
prob(item1, theta_n, derivative = 2)

theta <- rnorm(1)
item1 <- generate_item(model = "4PL")

# Probability of correct response
prob(item1, theta)

# First derivative of probability of correct response:
prob(item1, theta, derivative = 1)

# Second derivative of probability of correct response:
prob(item1, theta, derivative = 2)


# Multiple theta values
theta_n <- rnorm(5)

prob(item1, theta_n)
prob(item1, theta_n, derivative = 1)
prob(item1, theta_n, derivative = 2)

theta <- rnorm(1)
item1 <- generate_item(model = "GRM")

# Probability of correct response
prob(item1, theta)

# First derivative of probability of correct response:
prob(item1, theta, derivative = 1)

# Multiple theta values
theta_n <- rnorm(5)

prob(item1, theta_n)
prob(item1, theta_n, derivative = 1)


item4 <- generate_item(model = "GRM", n_categories = 5)
prob(item4, theta)

# Partial Credit Model
theta <- rnorm(1)
item1 <- generate_item(model = "PCM")

# Probability of correct response
prob(item1, theta)

# First derivative of probability of correct response:
prob(item1, theta, derivative = 1)

# Second derivative of probability of correct response:
prob(item1, theta, derivative = 2)


# Multiple theta values
theta_n <- rnorm(5)

prob(item1, theta_n)
prob(item1, theta_n, derivative = 1)
prob(item1, theta_n, derivative = 2)


item3 <- generate_item(model = "GPCM2", n_categories = 3)
prob(item3, theta)

theta <- rnorm(1)
item1 <- generate_item(model = "GPCM")

# Probability of correct response
prob(item1, theta)

# First derivative of probability of correct response:
prob(item1, theta, derivative = 1)

# Second derivative of probability of correct response:
prob(item1, theta, derivative = 2)


# Multiple theta values
theta_n <- rnorm(5)

prob(item1, theta_n)
prob(item1, theta_n, derivative = 1)
prob(item1, theta_n, derivative = 2)

# Probability of each response category for Generalized Partial Credit Model
item2 <- generate_item(model = "GPCM", n_categories = 4)
prob(item2, theta)

# First derivative of each response category
prob(item2, theta, derivative = 1)

# Second derivative of each response category
prob(item2, theta, derivative = 2)

theta <- rnorm(1)
item1 <- generate_item(model = "GPCM2")

# Probability of correct response
prob(item1, theta)

# First derivative of probability of correct response:
prob(item1, theta, derivative = 1)

# Second derivative of probability of correct response:
prob(item1, theta, derivative = 2)


# Multiple theta values
theta_n <- rnorm(5)

prob(item1, theta_n)
prob(item1, theta_n, derivative = 1)
prob(item1, theta_n, derivative = 2)


theta <- rnorm(1)
ip <- generate_ip(model = "3PL")

# Probability of correct response
prob(ip, theta)

# First derivative of probability of correct response:
prob(ip, theta, derivative = 1)

# Second derivative of probability of correct response:
prob(ip, theta, derivative = 2)

# Multiple theta
theta_n <- rnorm(3)
prob(ip, theta_n)
prob(ip, theta_n, derivative = 1)
prob(ip, theta_n, derivative = 2)

# Extract probabilities of correct response (i.e. response is "1")
sapply(prob(ip, theta_n), `[`, TRUE, "1")
# Probabilities of incorrect response
sapply(prob(ip, theta_n), `[`, TRUE, "0")

# Probability of each response category for Generalized Partial Credit Model
ip <- generate_ip(model = "GPCM", n = 4, n_categories = c(3, 4, 6, 5))
prob(ip, theta)

# First derivative of each response category
prob(ip, theta, derivative = 1)

# Second derivative of each response category
prob(ip, theta, derivative = 2)

# Probability of a mixture of items models
ip <- generate_ip(model = c("GPCM", "2PL", "3PL", "GPCM"),
                  n_categories = c(4, 2, 2, 3))
prob(ip, theta)

# Multiple theta
prob(ip, theta_n)

# Extract probabilities of score "2" for each theta value
sapply(prob(ip, theta_n), `[`, TRUE, "2")

theta <- rnorm(1)
t1 <- generate_testlet(model_items = "3PL")

# Probability of correct response
prob(t1, theta)

# First derivative of probability of correct response:
prob(t1, theta, derivative = 1)

# Second derivative of probability of correct response:
prob(t1, theta, derivative = 2)

}
\author{
Emre Gonulates
}
