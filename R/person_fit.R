
###############################################################################@
############################# person_fit (generic) #############################
###############################################################################@
#' Calculate person-fit indices
#' @description
#' \code{person_fit} calculates the fit of a person to a given psychometric
#' model.
#'
#'
#' @param ip An \code{\link{Item-class}}, \code{\link{Itempool-class}} or a
#'   \code{\link{Testlet-class}} object.
#' @param resp A vector of item responses.
#' @param theta An vector containing ability parameters.
#' @param type The type of the person-fit index.
#'
#' @return A vector of person-fit index values.
#'
#' @include item-class.R
#' @include itempool-class.R
#' @include item-class-methods.R
#' @include itempool-class-methods.R
#' @include response_set-class.R
#'
#' @author Emre Gonulates
#'
#' @rdname person_fit
#'
setGeneric("person_fit", function(resp, ip, theta, type = "lz")
{standardGeneric("person_fit")})


# ###############################################################################@
# ############################# person_fit (Item) ################################
# ###############################################################################@
# # @export
# #
# # @rdname person_fit
# #
# setMethod(
#   f = "person_fit", signature = c(ip = "Item"),
#   function(ip, resp, theta, type = "lz"){
#     if (type == "lz") {
#       # If the Model is one of the following: "irt1PM" "irt2PM" "irt3PM" "irt4PM"
#       if (ip$model %in% UNIDIM_DICHO_MODELS) {
#         # Calculate l_0:
#         l0 <- resp_loglik(ip = ip, resp = resp, theta = theta)
#         # Expected value of L_0
#         p <- prob(ip = ip, theta = theta)
#         q <- 1 - p
#         el0 <- p * log(p) + q * log(q)
#         # Variance of L_0
#         vl0 <- p * q * log(p/q)^2
#         return((l0 - el0) / sqrt(vl0))
#       } else
#         stop("Currently this method is only available for dichotomous IRT ",
#              "models.")
#     } else
#       stop("This method has not been implemented yet.")
#   }
# )


###############################################################################@
############################# person_fit (Response_set) ########################
###############################################################################@
#' @export
#'
#' @rdname person_fit
#'
setMethod(
  f = "person_fit", signature = c(resp = "Response_set", ip = "Itempool"),
  function(resp, ip, theta, type = "lz"){
    if (type == "lz") {
      return(lz_response_set_cpp(resp, theta, ip))
    } else
      stop("This method has not been implemented yet.")
  }
)

###############################################################################@
############################# person_fit (Itempool) ############################
###############################################################################@
#' @export
#'
#' @rdname person_fit
#'
setMethod(
  f = "person_fit", signature = c(ip = "Itempool"),
  function(resp, ip, theta, type = "lz"){
    if (type == "lz") {
      resp <- response_set(resp, ip = ip)
      return(lz_response_set_cpp(resp, theta, ip))
      # # If the Model is one of the following: "irt1PM" "irt2PM" "irt3PM" "irt4PM"
      # if (all(ip$model %in% UNIDIM_DICHO_MODELS)) {
      #   # Calculate l_0:
      #   l0 <- sum(resp_loglik(ip = ip, resp = resp, theta = theta),
      #             na.rm = TRUE)
      #   # Expected value of L_0
      #   p <- prob(ip = ip, theta = theta)
      #   q <- 1 - p
      #   el0 <- sum(p * log(p) + q * log(q), na.rm = TRUE)
      #   # Variance of L_0
      #   vl0 <- sum(p * q * log(p/q)^2, na.rm = TRUE)
      #   return((l0 - el0) / sqrt(vl0))
      # } else
      #   stop("Currently this method is only available for dichotomous IRT ",
      #        "models.")
    } else
      stop("This method has not been implemented yet.")
  }
)


###############################################################################@
############################# person_fit (Testlet) #############################
###############################################################################@
#' @export
#'
#' @rdname person_fit
#'
setMethod(
  f = "person_fit", signature = c(ip = "Testlet"),
  function(resp, ip, theta, type = "lz"){
    return(person_fit(ip = ip@item_list, resp = response_set(resp),
                      theta = theta, type = type))
  }
)
