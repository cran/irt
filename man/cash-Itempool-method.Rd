% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Itempool-class-methods.R
\name{$,Itempool-method}
\alias{$,Itempool-method}
\title{Get slots of the an \code{\link{Item-class}} object.}
\usage{
\S4method{$}{Itempool}(x, name)
}
\arguments{
\item{x}{An \code{\link{Itempool-class}} object from which to extract
element(s) or in which to replace element(s).}

\item{name}{Name of the parameter.
Available values:
\describe{
  \item{\strong{\code{'id'}}}{Extract \code{id}'s of all items
    and testlets.
    This will not extract the \code{id}'s of items within the
    testlet.}
  \item{\strong{\code{'content'}}}{Extract \code{content}'s of
    all items and testlets.
    This will not extract the \code{content}'s of items within the
    testlet.}
  \item{\strong{\code{'model'}}}{Extract \code{model}'s of
    all items and testlets.
    This will not extract the \code{model}'s of items within the
    testlet.}
  \item{\strong{\code{'misc'}}}{Extract \code{misc} parameters of
    all items and testlets.
    This will not extract the \code{misc} parameters of items
    within the testlet.}
  \item{\strong{\code{'item_list'}}}{Extract individual elements of
    item pool. If there are testlets in the item pool, a testlet
    will be an item of the resulting list. If individual items
    within the testlet is desired to be elements of the list, then
    use \code{$items}.}
  \item{\strong{\code{'items'}}}{Extract individual items
    within the item pool. If there are testlets in the item pool
    individual elements of the testlet will be extracted. Resulting
    list will only consist of \code{\link{Item-class}} objects.
  }
  \item{\strong{\code{'parameters'}}}{Extract \code{parameters}'s of
    all items and testlets.
    This will not extract the \code{parameters}'s of items within
    the testlet.}
  \item{\strong{\code{'se_parameters'}}}{Extract
    \code{se_parameters}'s of all items and testlets.
    This will not extract the \code{se_parameters}'s of items
    within the testlet.}
  \item{\strong{\code{'n'}}}{Return a list with three objects:
    \code{elements} the number of standalone items and testlets.
    \code{testlets} the number of Testlet objects.
    \code{items} the sum of the number of items within testlets and
    standalone items.
    }
  \item{\strong{\code{'max_score'}}}{Returns the maximum possible
    raw score of the item pool.
    }
  \item{\strong{\code{'resp_id'}}}{Extract
    \code{id}'s of all standalone items and items within the
    testlets. It will not return testlet \code{id}'s. This is
    mainly to get the \code{id}'s of items which has a response.
    }
  \item{\strong{\code{'resp_content'}}}{Extract
    \code{content}'s of all standalone items and items within the
    testlets. It will not return testlet \code{content}'s. This
    is mainly to get the \code{content}'s of items which has a
    response.
    }
  \item{\strong{\code{'resp_model'}}}{Extract
    \code{model}'s of all standalone items and items within the
    testlets. It will not return testlet \code{model}'s. This is
    mainly to get the \code{model}'s of items which has a response.
    }
  \item{\strong{\code{'resp_misc'}}}{Extract
    \code{misc} fields of all standalone items and items within
    the testlets. It will not return testlet \code{misc} fields.
    }
  \item{\strong{\code{'resp_item_list'}}}{Combine items that are
    not in a testlet and items within a testlet and return a list
    object. This list does not contain any \code{Testlet} objects.
    All of the elements are \code{Item} objects. If there are no
    testlets in the item pool, then this argument will be the
    same as \code{$item_list}.
    }
  \item{\strong{\code{'resp_max_score'}}}{Extract the maximum score
    each standalone item can get.
    }
}}
}
\value{
This operation will return a numeric object.
}
\description{
Get slots of the an \code{\link{Item-class}} object.
}
\examples{
item1 <- methods::new("Item", model =  "3PL", id = 'item23',
                      content = 'Geometry',
                      parameters = list(b = 2, c = .12, a = 1.2, D = 1))
item1$a
item1$D
item1$model
item1$id
item1$content
}
\author{
Emre Gonulates
}
