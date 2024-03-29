% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R
\name{classification_indices}
\alias{classification_indices}
\title{Calculate classification accuracy and consistency}
\usage{
classification_indices(
  method = "recursive",
  ip = NULL,
  theta = NULL,
  theta_cs = NULL,
  raw_cs = NULL,
  resp = NULL,
  se = NULL,
  perf_categories = NULL,
  n_theta = 100,
  theta_lower_bound = -6,
  theta_upper_bound = 6,
  cat_labels = NULL
)
}
\arguments{
\item{method}{The method of classification accuracy and consistency
 calculation method. Following methods are available:

  \describe{
    \item{\code{'rudner'}}{
      Rudner (2000, 2005) based classification accuracy and consistency
      indices.

      Following values should be provided for this method:
      \code{theta, se, theta_cs}.
      Following values can optionally be provided for this method:
      \code{perf_categories, cat_labels}.
    }
    \item{\code{'guo'}}{Guo (2006) based classification accuracy and
      consistency indices.

      Note that calculation times can be long for this method. The value of
      \code{n_theta} can be decreased to speed up the function but this
      will reduce the accuracy of the index.

      Following values should be provided for this method:
      \code{ip, resp, theta_cs} and either one of \code{theta} or
      \code{perf_categories}.
      Following values can optionally be provided for this method:
      \code{n_theta, theta_lower_bound, theta_upper_bound, cat_labels}.
    }
    \item{\code{'recursive'}}{Lee (2010) based classification
      accuracy and consistency indices.

      Following values should be provided for this method:
      \code{ip, theta} and either one of these \code{theta_cs, raw_cs}.
      Following values can optionally be provided for this method:
      \code{perf_categories, cat_labels}.
    }
  }}

\item{ip}{An \code{\link{Itempool-class}} object. Item pool parameters can be
composed of any combination of unidimensional dichotomous or polytomous
items. Required for \code{"guo"} and \code{"recursive"} methods.}

\item{theta}{A numeric vector representing the abilities of examinees.
Required for \code{'rudner'} and \code{'recursive'} method. For
\code{"guo"} method, this vector will be used to get performance category
of each examinee if \code{perf_categories} is \code{NULL}. The default
value is \code{NULL}. For \code{method = "guo"} either \code{theta} or
\code{perf_categories} should be provided.}

\item{theta_cs}{A sorted (ascending order) numeric vector representing the
theta scale cut scores. Do not include \code{-Inf} or \code{Inf}. Required
for \code{'rudner'} and \code{'guo'} method; required for 'recursive' if
\code{raw_cs} is not provided.}

\item{raw_cs}{A sorted (ascending order) numeric vector of summed-score cut
score values. Do not include 0 or the maximum possible score of the test in
this vector. Required for 'recursive' method if 'theta_cs' is not provided.}

\item{resp}{A \code{\link{Response_set-class}}, a \code{matrix} or a
\code{data.frame} object that holds responses. If \code{matrix} or a
\code{data.frame} provided, they will be converted to a
\code{\link{Response_set-class}}. Required for \code{'guo'} method.}

\item{se}{A numeric vector representing the standard errors of ability
estimates. Required for \code{'rudner'} method.}

\item{perf_categories}{An integer  vector representing the performance
categories of examinees. The number 1 should represent the lowest category.
For example if there are three cut scores the valid values can only be: 0,
1, 2 and 3. This  vector will be used \code{theta} is \code{NULL}. The
default value is \code{NULL}. Either \code{theta} or \code{perf_categories}
should be provided. Can optional be provided for all methods.}

\item{n_theta}{An integer representing the number of equally spaced theta
points between cut scores. The default value is 100. Use larger values to
increase accuracy but larger numbers will also slow the speed of
calculation. Can optionally be provided for the \code{'guo'} method.}

\item{theta_lower_bound}{A number representing the lower bound for cut
scores. The default value is -6. Can optionally be provided for the
\code{'guo'} method.}

\item{theta_upper_bound}{A number representing the upper bound for cut
scores. The default value is 6. Can optionally be provided for the
\code{'guo'} method.}

\item{cat_labels}{A string vector representing the labels of the categories.
The length of the vector should be one more than the length of the cut
scores. The default value is \code{NULL} where the categories will be
labeled as 1, 2, ..., (number of cut scores plus one). For example, if
there are three cut scores category labels can be:
\code{c("Unsatisfactory", "Basic", "Mastery", "Advanced")}. Can optional be
provided for all methods.}
}
\value{
A list of following elements:
  \describe{
    \item{\code{category_prob}}{A numeric vector representing the
      performance category classification probabilities of each examinee.}
    \item{\code{ca}}{Marginal (overall) classification accuracy index}
    \item{\code{cc}}{Marginal (overall) classification consistency index}
    \item{\code{ind_cs_ca}}{Individual cut score classification accuracy
      indices. This value will only be calculated when there are more than
      one cut score.}
    \item{\code{ind_cs_ca}}{Individual cut score classification consistency
      indices. This value will only be calculated when there are more than
      one cut score.}
  }
}
\description{
Calculate classification accuracy and consistency
}
\examples{

ip <- generate_ip(model = sample(c("GPCM", "2PL"), 20, TRUE))
n_examinee <- 100

true_theta <- rnorm(n_examinee)
resp_set <- generate_resp_set(ip = ip, theta = true_theta, prop_missing = .2)
theta_est <- est_ability(resp = resp_set, ip = ip, method = "eap")
se <- theta_est$se
theta_est <- theta_est$est
raw_score <- est_ability(resp = resp_set, method = "sum_score")$est

# Cut score
theta_cs <- c(-1, 0, 1.5)
raw_cs <- round(rsss(ip = ip, scale_score = theta_cs))

# Rudner (2000, 2005) based indices:
classification_indices(method = "rudner", theta = theta_est, se = se,
                       theta_cs = theta_cs)

# Guo (2006) based indices:
classification_indices(method = "guo", ip = ip, resp = resp_set,
                       theta = theta_est, theta_cs = theta_cs)

# Recursive method based indices:
classification_indices(method = "recursive", ip = ip, theta = theta_est,
                       theta_cs = theta_cs)
# Use raw score cut scores with recursive method
classification_indices(method = "recursive", ip = ip, theta = theta_est,
                       raw_cs = raw_cs)

}
\references{
Guo, F. (2006). Expected classification accuracy using the latent
distribution. Practical Assessment, Research, and Evaluation, 11(1), 6.

Lee, W. C. (2010). Classification consistency and accuracy for complex
assessments using item response theory. Journal of Educational Measurement,
47(1), 1-17.

Rudner, L. M. (2000). Computing the expected proportions of misclassified
examinees. Practical Assessment, Research, and Evaluation, 7(1), 14.

Rudner, L. M. (2005). Expected classification accuracy. Practical Assessment,
Research, and Evaluation, 10(1), 13.

Wyse, A. E., & Hao, S. (2012). An evaluation of item response theory
classification accuracy and consistency indices. Applied Psychological
Measurement, 36(7), 602-624.
}
\author{
Emre Gonulates
}
