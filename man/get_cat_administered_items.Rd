% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cat_sim_helper_functions.R
\name{get_cat_administered_items}
\alias{get_cat_administered_items}
\title{Get administered items from a CAT output}
\usage{
get_cat_administered_items(cat_sim_output)
}
\arguments{
\item{cat_sim_output}{This is a list object containing elements that are
"cat_output" class.}
}
\value{
For \code{cat_output} with only one adaptive test, an
  \code{Itempool} class object will be returned. For \code{cat_output} with
  more than one adaptive tests, a list of \code{Itempool} class objects will
  be returned.
}
\description{
This function returns an item pool object of the
  administered items using the items in estimate history. If there is one
}
\examples{
cd <- create_cat_design(ip = generate_ip(n = 30), next_item_rule = 'mfi',
                        termination_rule = 'max_item',
                        termination_par = list(max_item = 10))
cat_data <- cat_sim(true_ability = rnorm(10), cd = cd)
get_cat_administered_items(cat_data)
}
\author{
Emre Gonulates
}
