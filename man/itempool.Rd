% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/itempool-class-methods.R
\name{itempool}
\alias{itempool}
\title{Create an \code{Itempool} object}
\usage{
itempool(...)
}
\arguments{
\item{...}{The object that is desired to be converted to an 'Itempool'
object. Also additional arguments related to the Itempool.}
}
\value{
An \code{\link{Itempool-class}} object.
}
\description{
This method creates a new \code{\link{Itempool-class}} object.
}
\examples{
# Create an item pool with two 2PL items
itempool(a = c(1, 1.4), b = c(-2, 1))
itempool(a = c(1, 1.4), b = c(-2, 1), model = "2PL")
# Set D parameter
itempool(a = c(1, 1.4), b = c(-2, 1), D = 1.7)
# Set item IDs
itempool(a = c(1, 1.4), b = c(-2, 1), item_id = c("i1", "i2"))
# Set content
itempool(a = c(1, 1.4), b = c(-2, 1), content = c("Algebra", "Geometry"))

# Create 3PL items from a data frame:
ipdf <- data.frame(a = c(.9, .79, 1.26),
                   b = c(-1, .43, -2.3),
                   c = c(.2, .38, .25))
itempool(ipdf)

# Create GRM (Graded Response Model) items from a data frame
ipdf <- data.frame(a = rlnorm(10, 0, .3), b1 = rnorm(10), b2 = rnorm(10))
itempool(ipdf, model = "GRM")

# Create a Rasch model item pool
itempool(b = c(-1, 0.2, 1.1), model = "Rasch")

# Add 'misc' field:
ip <- itempool(b = rnorm(2), item_id = paste0("t1-i", 1:2),
               misc = list(list(sympson_hetter_k = .8),
                           list(sympson_hetter_k = .9)))
ip[[1]]  # First item of the item pool

}
\author{
Emre Gonulates
}
