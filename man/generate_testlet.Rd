% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_objects.R
\name{generate_testlet}
\alias{generate_testlet}
\title{Generate a random Testlet object}
\usage{
generate_testlet(
  model = "BTM",
  n = NULL,
  item_models = "3PL",
  item_id_preamble = NULL,
  n_categories = 4,
  ...
)
}
\arguments{
\item{model}{The model of the Testlet}

\item{n}{The number of items in the Testlet.}

\item{item_models}{A single model name or a vector of model names with the
size of n that represents the models of items in the Testlet object.}

\item{item_id_preamble}{The preamble for the item ids within the Testlet.}

\item{n_categories}{For polytomous items, designate the number of categories
each item should have. It can be a single integer value larger than 1. In
this case all of the polytomous items of the testlet will have this number
of categories. It can be a vector of length \code{n} designating the
categories of each item. For dichotomous items, the values in
\code{n_categories} will be ignored.}

\item{...}{Additional parameters passed to \code{testlet()} function.}
}
\value{
A \code{\link{Testlet-class}} object
}
\description{
Generate a random Testlet object
}
\examples{
# By default, a Testlet object with '3PL' model items generated
generate_testlet()
# Designate the number of items in the testlet
generate_testlet(n = 12)
# Set the ID of the testlet
generate_testlet(testlet_ = "my-testlet")
# Designate the ID of testlet and preamble for item ids
generate_testlet(testlet_id = "my-testlet", item_id_preamble = "mt-")
# Generate item pools for other models
generate_testlet(item_model = "Rasch")
generate_testlet(item_model = "1PL")
generate_testlet(item_model = "2PL")
generate_testlet(item_model = "4PL")
generate_testlet(item_model = "GRM") # Graded Response Model
generate_testlet(item_model = "GPCM") # Generalized Partial Credit Model
generate_testlet(item_model = "PCM") # Partial Credit Model
generate_testlet(item_model = "GPCM2") # Reparameterized GPCM
# Mixture of models
generate_testlet(item_models = c("4PL", "Rasch"))
generate_testlet(model = c("2PL", "GRM", "Rasch"), n = 11)

# Generating multiple testlet objects with custom ids
sapply(paste0("testlet-", 1:4), function(x) generate_testlet(testlet_id = x))

# Generate testlet with dichotomous and polytomous with different number of
# categories.
generate_testlet(
  item_models = c("3PL", "GRM", "GPCM", "GRM", "2PL"),
  n_categories = c(2, 3, 6, 7, 2))

# # Generating multiple testlet objects with custom ids and item models and
# # put them in an item pool:
# temp_list <- list(ids = paste0("testlet-", 1:3),
#                   item_models = c("Rasch", "2PL", "GPCM"))
# itempool(sapply(1:length(temp_list$item_id), function(i)
#   generate_testlet(item_id = temp_list$item_id[i],
#   item_models = temp_list$item_models[i])))

}
\author{
Emre Gonulates
}
