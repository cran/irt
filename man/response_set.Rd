% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/response_set-class-methods.R
\name{response_set}
\alias{response_set}
\title{Create \code{\link{Response_set-class}} object}
\usage{
response_set(
  x,
  data_format = "wide",
  ip = NULL,
  examinee_id_var = NULL,
  testlet_id_var = NULL,
  item_id_var = NULL,
  score_var = NULL,
  raw_response_var = NULL,
  order_var = NULL,
  response_time_var = NULL,
  misc_var = NULL,
  misc_unique_var = NULL,
  misc = NULL,
  fill_na_score = NULL
)
}
\arguments{
\item{x}{A \code{matrix} or \code{data.frame} holding item scores. See the
description about the options. Additionally, it can be a list of
\code{\link{Response-class}} objects.}

\item{data_format}{A string value representing the format of the data
  \code{x} supplied. The default value is \code{"wide"}.
  The following options are available:
  \describe{
    \item{"wide"}{\code{x} can be in wide format data where a \code{matrix}
      or \code{data.frame} where rows represents examinees and columns
      represent items. Each row will be converted to a
      \code{\link{Response-class}} object.

      If the columns has names (and an \code{\link{Itempool-class}} object
      has not been supplied), then the \code{item_id}s will be supplied by
      the column names. If neither column names nor an
      \code{\link{Itempool-class}} object supplied, default \code{item_id}s
      will be given.

      If rows has names, those will be used as \code{examinee_id}s.
    }
    \item{"long"}{\code{x} can be in long format where \code{data.frame}
      with at least three columns: (1) a column for \code{examinee_id}, (2) a
      column for \code{item_id} and (3) a column for either \code{score}s or
      \code{raw_response}s. Additional columns can be added such as
      \code{testlet_id}, item \code{order}, \code{response_time}.
    }
  }}

\item{ip}{Optionally an \code{\link{Itempool-class}} object that is holding
the item parameters can be supplied to check whether Response_set object
created is compatible with the \code{\link{Itempool-class}} object.}

\item{examinee_id_var}{A string for the column name that holds examinee
ids, if \code{x} is in long format.}

\item{testlet_id_var}{A string for the column name that holds testlet
ids, if \code{x} is in long format.}

\item{item_id_var}{A string for the column name that holds item
ids, if \code{x} is in long format.}

\item{score_var}{A string for the column name that holds examinee
scores, if \code{x} is in long format.}

\item{raw_response_var}{A string for the column name that holds raw
responses of the examinees, if \code{x} is in long format.}

\item{order_var}{A string for the column name that holds the administration
order of items, if \code{x} is in long format.}

\item{response_time_var}{A string for the column name that holds response
time information of the items, if \code{x} is in long format.}

\item{misc_var}{A string for the column names that are holding the
miscellaneous information of the items. Available only when \code{x} is
in long format. Within an examinee, if there is additional information
for each item (for example, item's type, item's reading level, examinee's
raw response to an item, whether an item is operational or not,
the date/time item is administered, ratings of multiple raters, etc.),
in the dataset, this information can be passed. Later in the code,
such information can be extracted by \code{$} operator. See examples.}

\item{misc_unique_var}{A string for the column names that are holding the
miscellaneous information of the items. Different than \code{misc_var},
these columns are assumed to be the same within an examinee, so only the
unique value of this column within an examinee will be saved. Examples of
variables for this column is gender, race, ability score, school of the
examinee that will not vary from one item to another within an examinee.
The argument is only available when \code{data_format = "long"}.}

\item{misc}{A list of miscellaneous variables that needs to be added to
the \code{Response_set} object.}

\item{fill_na_score}{If some examinees do not answer all items, the value
  \code{fill_na_score} will be replaced by the scores of unanswered items.
  If an \code{ip} value provided, 'all items' will be all of the items in
  the item pool. Otherwise, all items will be the list of all unique
  \code{item_id} values.

  Currently, this feature only works when \code{x} is a data frame or matrix.}
}
\value{
A \code{\link{Response_set-class}} object.
}
\description{
This function creates a \code{\link{Response_set-class}} object
  from various types of data sets. Currently following scenarios are
  supported:
}
\examples{
##### Wide format data #####
## Example 1
x_wide <- matrix(sample(0:1, 35, TRUE), nrow = 7, ncol = 5)
response_set(x_wide)

## Example 2
ip <- generate_ip(n = 6)
# simulate responses for 10 examinees
resp_matrix <- sim_resp(ip = ip, theta = rnorm(10), prop_missing = .2,
                        output = "matrix")
# convert it to tibble
resp_wide <- as.data.frame(resp_matrix)
resp_wide$stu_id <- rownames(resp_matrix)
# Create a Response_set object:
resp_set <- response_set(resp_wide, data_format = "wide", ip = ip,
                         examinee_id_var = "stu_id")
# Retrieve examinee ids:
resp_set$examinee_id
# Fourth examinee:
resp_set[[4]]
# Scores of 6th examinee
resp_set[[6]]$score


##### Long format data #####
x_long <- data.frame(examinee_id = c("stu1", "stu1", "stu1", "stu2", "stu2"),
                     item_id = c("i1", "i2", "i4", "i1", "i2"),
                     scr = c(0, 1, 0, 1, 0),
                     rwscore = c("A", "D", "B", "C", "D"),
                     resptime = c(33, 55, 22, 66, 31),
                     # These will be passed to misc
                     item_type = c("MC", "MC", "MS", "SA", "MC"),
                     lexile_level = c(1, 4, 3, 2, 1),
                     word_count = c(123, 442, 552, 342, 666),
                     ability = c(1.1, 1.1, 1.1, -.2, -.2),
                     grade = c("7", "7", "7", "11", "11")
                     )

resp_set <- response_set(x = x_long,
                         data_format = "long",
                         examinee_id_var = "examinee_id",
                         item_id_var = "item_id",
                         score_var = "scr",
                         raw_response_var = "rwscore",
                         response_time_var ="resptime",
                         misc_var = c("item_type", "lexile_level"),
                         misc_unique_var = c("ability", "grade")
                         )

resp_set[[1]]  # Response of the first examinee
resp_set$item_type  # extract item_type of each examinee
resp_set$grade  # extract grade of each examinee

# Also, additional examinee level miscellaneous information can be added:
resp_set$gender <- c("M", "F")
resp_set[[2]]$gender  # access second examinee's gender.
resp_set$gender

# Fill missing values with 0.
response_set(x = x_long,
             data_format = "long",
             examinee_id_var = "examinee_id",
             item_id_var = "item_id",
             score_var = "scr",
             raw_response_var = "rwscore",
             response_time_var ="resptime",
             misc_var = c("item_type", "lexile_level"),
             fill_na_score = 0
             )
}
\author{
Emre Gonulates
}
