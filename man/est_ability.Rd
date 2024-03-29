% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ability_estimation.R
\name{est_ability}
\alias{est_ability}
\title{Estimate Examinee Ability}
\usage{
est_ability(
  resp,
  ip = NULL,
  method = c("eap", "ml", "map", "bm", "owen", "sum_score"),
  ...,
  prior_dist = c("norm", "unif", "lnorm", "gamma", "t", "cauchy"),
  prior_pars = c(0, 1),
  theta_range = c(-5, 5),
  number_of_quads = 41,
  tol = 1e-06,
  output_type = c("list", "data.frame", "tibble")
)
}
\arguments{
\item{resp}{A \code{\link{Response_set-class}}, \code{matrix} or a
\code{data.frame} object holding responses. Missing responses are excluded
from the ability estimation.}

\item{ip}{An \code{\link{Item-class}}, \code{\link{Itempool-class}} or a
\code{\link{Testlet-class}} object. If \code{ip} is not an
\code{\link{Itempool-class}} object, the function attempts to convert it.
While the default is \code{NULL}, this argument is required for all methods
except when \code{method = "sum_score"}.}

\item{method}{The method used for ability estimation. The default is
  \code{"eap"}.

  Available methods:
  \describe{
    \item{\strong{\code{'sum_score'}}}{Basic sum (raw) score of responses.}
    \item{\strong{\code{'owen'}}}{Owen's Bayesian Ability Estimation.

      This method is suitable for dichotomous IRT models (e.g., 'Rasch',
      '1PL', '2PL', '3PL' and '4PL'). Testlet groupings are ignored and items
      within testlets are treated as standalone items.

      Formulas were implemented in Owen (1975) and Vale (1977). The original
      formulation does not include the D parameter. If \code{D = 1}, the
      original solution is obtained. If \code{D = 1.7}, the \code{a}
      parameter is multiplied by this number.

      The user needs to provide prior parameters, i.e., \code{prior_pars}.
      These should be a numeric vector of length two, with the first
      component as the prior mean and the second as the prior standard
      deviation (not variance). For example, if the prior mean is 0.1 and the
      prior standard deviation is 2, set the prior parameters as
      \code{prior_pars = c(0.1, 2)}.
      }

    \item{\strong{\code{'ml'}}}{Maximum Likelihood Ability Estimation
      via Newton-Raphson Algorithm.}
    \item{\strong{\code{'eap'}}}{Expected-a-Posteriori Ability Estimation.
      Prior information must be provided for this function. The number of
      quadrature points can also be specified using the argument
      \code{number_of_quads}. }
    \item{\strong{\code{'map'}} or \strong{\code{'bm'}}}{Maximum-a-Posteriori
      Ability Estimation (or Bayes Modal estimation). Prior information must
      be provided for this function. Currently, only \code{'norm'} prior
      distribution is available.}
  }}

\item{...}{Additional arguments passed to specific methods.}

\item{prior_dist}{The shape of the prior distribution. Available options are:
\describe{
  \item{'norm'}{Normal distribution}
  \item{'unif'}{Uniform distribution}
  \item{'t'}{t distribution}
  \item{'cauchy'}{Cauchy distribution}
}
The default value is \code{'norm'}.}

\item{prior_pars}{Parameters of the prior distribution. Default value is
  \code{c(0, 1)}, where 0 is the mean and 1 is the standard deviation of the
  default normal prior distribution. For example, uniform prior parameter can
  be set as \code{c(a, b)} where \code{a} is the minimum value and \code{b}
  is the maximum value. For \code{t} distribution, prior parameter can be set
  as \code{df} to represent the degree of freedom. For Cauchy distribution,
  prior parameters can be set as \code{c(location, scale)}.

  If method is \code{"owen"}, provide \code{c(<Prior Mean>, <Prior SD>)}.}

\item{theta_range}{The limits of the ability estimation scale. The estimation
result will be bounded within this interval. Default is \code{c(-5, 5)}.}

\item{number_of_quads}{Number of quadratures. The default value is 41. As
this number increases, the precision of the estimate will also increase.}

\item{tol}{The precision level of ability estimate. The final ability
estimates will be rounded to remove precision smaller than the \code{tol}
value. Default is \code{1e-06}.}

\item{output_type}{A string specifying the output type of the function.
Default is \code{"list"}. Options include:
\describe{
  \item{"list"}{Function returns a \code{list} object with elements
    \code{est} and \code{se}.}
  \item{"data.frame"}{Function returns a \code{data.frame} object with
    columns \code{examinee_id}, \code{est} and \code{se}.}
  \item{"tibble"}{If the \code{tibble} package is available, the function
    returns a \code{tibble} object with columns \code{examinee_id},
    \code{est} and \code{se}.}
}}
}
\value{
\code{est} The estimated examinee abilities. If the response vector
  for a subject contains all \code{NA}s, then \code{est} will be \code{NA} to
  differentiate from cases where all answers are incorrect.

\code{se} The standard errors of the ability estimates. For
  \code{"sum_score"} method, all standard errors will be \code{NA}. For
  Bayesian methods (like EAP, MAP or Owen's), this value is the square root
  of the posterior variance.
}
\description{
This function estimates examinee ability using different methods, including
Owen's Bayesian estimation, Maximum Likelihood estimation,
Maximum-a-Posteriori and Expected-a-Posteriori.
}
\examples{
ip <- generate_ip(n = 7)
resp <- sim_resp(ip, theta = rnorm(3))

### EAP estimation ###
est_ability(resp, ip)
est_ability(resp, ip, number_of_quads = 81)
# The default prior_dist is 'norm'. prior_pars = c(mean, sd)
est_ability(resp, ip, prior_pars = c(0, 3))
# prior_pars = c(min, max)
est_ability(resp, ip, prior_dist = 'unif',  prior_pars = c(-3, 3))
# prior_pars = c(df)
est_ability(resp, ip, prior_dist = 't',  prior_pars = 3)
# prior_pars = c(location, scale)
est_ability(resp, ip, prior_dist = 'cauchy',  prior_pars = c(0, 1))


### MAP estimation (Bayes Modal estimation) ###
est_ability(resp, ip, method = "map")
# The default prior_dist is 'norm'. prior_pars = c(mean, sd)
est_ability(resp, ip, method = "map", prior_pars = c(0, 2))


### Maximum Likelihood estimation ###
est_ability(resp, ip, method = 'ml')
est_ability(resp, ip, method = 'ml', tol = 1e-8)
est_ability(resp = rep(1, length(ip)), ip, method = 'ml')
est_ability(resp = rep(1, length(ip)), ip, method = 'ml',
            theta_range = c(-3, 3))

### Owen's Bayesian ability estimation ###
est_ability(resp, ip, method = 'owen')
est_ability(resp, ip, method = 'owen', prior_pars = c(0, 3))



}
\references{
Owen, R. J. (1975). A Bayesian sequential procedure for quantal response in
the context of adaptive mental testing. Journal of the American Statistical
Association, 70(350), 351-356.

Vale, C. D., & Weiss, D. J. (1977). A Rapid Item-Search Procedure for
Bayesian Adaptive Testing. Research Report 77-4. Minneapolis, MN.
}
\author{
Emre Gonulates
}
