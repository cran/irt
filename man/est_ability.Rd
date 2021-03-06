% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ability_estimation.R
\name{est_ability}
\alias{est_ability}
\title{Ability Estimation of an examinee}
\usage{
est_ability(
  ip,
  resp,
  method = "eap",
  ...,
  prior_dist = "norm",
  prior_pars = c(0, 1),
  theta_range = c(-5, 5),
  number_of_quads = 41,
  tol = 1e-06
)
}
\arguments{
\item{ip}{An \code{\link{Item-class}}, \code{\link{Itempool-class}} or a
\code{\link{Testlet-class}} object.}

\item{resp}{A vector containing examinee responses. If there are missing
responses, they will not be included in the ability estimation.}

\item{method}{The method that will be used to estimate the ability.
  The default value is \code{"eap"}.

  Current methods are:
  \describe{
    \item{\strong{\code{'sum_score'}}}{Basic sum (raw) score of
      responses.}
    \item{\strong{\code{'owen'}}}{Owen's Bayesian Ability Estimation.

      This estimation method can be used only for dichotomous IRT
      models, 'Rasch', '1PL', '2PL', '3PL' and '4PL'.

      Formulas were implemented in Owen (1975) and Vale (1977).  Original
      formulation does not contain D parameter. If \code{D = 1} original
      solution will be obtained. If \code{D = 1.7}  the \code{a} parameter
      will be multiplied with this number.

      }
    \item{\strong{\code{'ml'}}}{Maximum Likelihood Ability Estimation
      via Newton-Raphson Algorithm}
    \item{\strong{\code{'eap'}}}{Expected-a-Posteriori Ability
      Estimation}
  }}

\item{...}{Additional arguments passed to specific methods}

\item{prior_dist}{The shape of the prior distribution. Currently following
distributions can be specified:
\describe{
  \item{'norm'}{Normal distribution}
  \item{'unif'}{Uniform distribution}
  \item{'t'}{t distribution}
  \item{'cauchy'}{Cauchy distribution}
}
Default value is \code{'norm'}.}

\item{prior_pars}{Parameters of the prior distribution. Default value is
  \code{c(0, 1)} where 0 is the mean and 1 is the standard deviation of the
  default prior distribution which is normal distribution. Also, for example,
  uniform prior parameter can be set as \code{c(a, b)} where \code{a} is the
  minimum value and \code{b} is the maximum value. For \code{t} distribution,
  prior parameter can be set as \code{df} to represent the degree of freedom.
  For Cauchy distribution, prior parameters can be set as \code{c(location,
  scale)}.

  If method is \code{"owen"}, provide \code{c(<Prior Mean>, <Prior SD>)}.}

\item{theta_range}{The limits of the ability estimation scale. The estimation
result will be limited to this interval. The default is \code{c(-5, 5)}.}

\item{number_of_quads}{Number of quadratures. The default value is 41. As
this number increases, the precision of the estimate will also increase.
The default value is \code{41}.}

\item{tol}{The precision level of ability estimate. The final ability
estimates will be rounded to remove the precision that is smaller than the
\code{tol} value. The default value is \code{1e-06}.}
}
\value{
\code{est} The ability estimated. If the response vector for a
  subject contains all \code{NA}s, then, in order to differentiate all
  incorrect and all NA, the \code{est} returned will be NA.

\code{se} The standard error(s) of the ability estimate(s). For
  \code{"sum_score"} method, all of the standard errors will be \code{NA}.
}
\description{
\code{est_ability} estimates ability using on various methods such as
Owen's Bayesian estimation, Maximum Likelihood estimation,
Expected-a-Posteriori.
}
\examples{
ip <- generate_ip()
resp <- sim_resp(ip, theta = rnorm(1))

# EAP estimation
est_ability(ip, resp)
est_ability(ip, resp, number_of_quads = 81)
est_ability(ip, resp, prior_pars = c(0, 3))
est_ability(ip, resp, prior_dist = 'unif',  prior_pars = c(-3, 3))
est_ability(ip, resp, prior_dist = 't',  prior_pars = 3)
est_ability(ip, resp, prior_dist = 'cauchy',  prior_pars = c(0, 1))

# Maximum Likelihood estimation
est_ability(ip, resp, method = 'ml')
est_ability(ip, resp, method = 'ml', tol = 1e-8)
est_ability(ip, resp = rep(1, length(ip)), method = 'ml')
est_ability(ip, resp = rep(1, length(ip)), method = 'ml',
            theta_range = c(-3, 3))

# Owen's Bayesian ability estimation
est_ability(ip, resp, method = 'owen')
est_ability(ip, resp, method = 'owen', prior_pars = c(0, 3))

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
