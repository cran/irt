#include <Rcpp.h>
#include "misc.h"
#include "itempool_class_methods.h"
#include "response_set_methods.h"
#include "resp_lik.h"
#include "resp_loglik.h"
using namespace Rcpp;


//#############################################################################@
//#############################################################################@
//########################### est_ability_4pm_nr ###############################
//#############################################################################@
//#############################################################################@

//#############################################################################@
//########################### est_ability_4pm_nr_itempool_single_iv_cpp ########
//#############################################################################@
double est_ability_4pm_nr_itempool_single_iv_cpp(
  Rcpp::NumericVector resp, Rcpp::S4 ip, Rcpp::NumericVector theta_range,
  double initial_est = 0, double criterion = 0.001) {
  double est = initial_est;
  double difference = criterion + 1;
  double firstDerivative, secondDerivative, adjustment, newEst;
  double minEst = theta_range[0];
  double maxEst = theta_range[1];

  if ((initial_est <= minEst) || (initial_est >= maxEst)) {
    est = 0;
  }
  while (difference > criterion) {
    firstDerivative = resp_loglik_bare_itempool_cpp(resp, est, ip, 1);
    secondDerivative = resp_loglik_bare_itempool_cpp(resp, est, ip, 2);

    // fabs = std::abs() function that calculates the absolute value.
    // -fabs(secondDerivative) is an adjustment to make algorithm to
    // the move correct direction
    adjustment = firstDerivative / (-fabs(secondDerivative));
    // Limit the move size by .5.
    if (fabs(adjustment) < .5) {
      newEst = est - adjustment;
    } else {
      newEst = est - 0.5 * ((adjustment > 0) - (adjustment < 0));
    }
    difference = fabs(newEst - est);
    // If estimate is out of bounds it will bound it.
    if ((newEst <= minEst) || (newEst >= maxEst)) {
      if ((est <= minEst) || (est >= maxEst))
      {
        if (est <= minEst) {
          return(minEst);
        } else {
          return(maxEst);
        }
      }
    }
    est = newEst;
  }
  return est;
}

//#############################################################################@
//########################### est_ability_4pm_nr_itempool_cpp ##################
//#############################################################################@
// [[Rcpp::export]]
double est_ability_4pm_nr_itempool_cpp(
  Rcpp::NumericVector resp, Rcpp::S4 ip, Rcpp::NumericVector theta_range,
  double criterion = 0.001,
  Rcpp::Nullable<Rcpp::NumericVector> initial_estimates = R_NilValue) {
  // Rcout << "est_ability_ml_nr - Stage 1" << std::endl;
  Rcpp::NumericVector init_est(3); // vector holding initial values
  double est; // Final estimate value that will be returned
  // Vector that holds estimates obtained using different initial values.
  bool all_equal = true;

  // Set the initial values if a NULL value provided
  if (initial_estimates.isNotNull()) {
    init_est = as<Rcpp::NumericVector>(initial_estimates);
  } else {
    init_est[0] = theta_range[0] + 2 * criterion;
    init_est[1] = 0;
    init_est[2] = theta_range[1] - 2 * criterion;
    // init_est = {theta_range[0] + 2 * criterion, 0,
    //             theta_range[1] - 2 * criterion};
  }
  int n_init_values = init_est.size();
  // Vector that holds estimates obtained using different initial values.
  Rcpp::NumericVector estimates(n_init_values);

  // Rcout << "est_ability_ml_nr - Stage 2" << std::endl;

  // Make sure that init_est has at least two values:
  if (init_est.size() < 2)
    stop("Please proivde at least two different initial values.");
  // Find ML ability estimates based on different initial values.

  for(int i = 0; i < n_init_values; i++) {
    estimates[i] = est_ability_4pm_nr_itempool_single_iv_cpp(
      resp, ip, theta_range, init_est[i], criterion);
    if (i > 0 && (std::fabs(estimates[i] - estimates[i-1]) >
                    n_init_values*criterion))
      all_equal = false;
    // Rcout << i << " - " << estimates[i] << std::endl;
  }
  est = estimates[0];
  // Rcout << "est_ability_ml_nr - Stage 5 - " << all_equal << " ; est = " <<
  // est << std::endl;

  // Check whether all elements of the vector are equal, if not return the
  // estimate value with the highest response log-likelihood
  if (!all_equal) {
    //Rcout << "est_ability_ml_nr - Stage 5.1 - Difference found!" << std::endl;

    // Calculate response log likelihood value at the estimate that was
    // derived from first initial estimate.
    double resp_ll = 0;
    double resp_ll_max = resp_loglik_bare_itempool_cpp(resp, estimates[0], ip,
                                                       0);

    for(int i = 0; i < n_init_values; i++) {
      // Rcout << "est_ability_ml_nr - Stage 5.2" << std::endl;
      resp_ll = resp_loglik_bare_itempool_cpp(resp, estimates[i], ip, 0);
      if (resp_ll > resp_ll_max) {
        est = estimates[i];
        resp_ll_max = resp_ll;
      }
    }
  }
  // Rcout << "est_ability_ml_nr - Stage 7 - " << est << std::endl;
  return(est);
}

//#############################################################################@
//########################### est_ability_4pm_nr_iv_response_cpp ###############
//#############################################################################@
//' This function estimates ability using Newton-Raphson method for a Response
//' object.
//' @noRd
//'
double est_ability_4pm_nr_iv_response_cpp(
    Rcpp::S4 resp,
    Rcpp::S4 ip,
    Rcpp::NumericVector theta_range = NumericVector::create(-5, 5),
    double initial_est = 0,
    double criterion = 0.001) {
  double est = initial_est;
  double difference = criterion + 1;
  double firstDerivative, secondDerivative, adjustment, newEst;
  double minEst = theta_range[0];
  double maxEst = theta_range[1];

  // Rcout << "  est_ability_ml_nr_iv - Stage 1" << std::endl;

  if ((initial_est <= minEst) || (initial_est >= maxEst)) {
    est = 0;
  }
  while (difference > criterion) {

    // Rcout << "  est_ability_ml_nr_iv - Stage 1.1" << std::endl;

    firstDerivative = resp_loglik_response_cpp(est, resp, ip, 1);
    secondDerivative = resp_loglik_response_cpp(est, resp, ip, 2);

    // fabs = std::abs() function that calculates the absolute value.
    // -fabs(secondDerivative) is an adjustment to make algorithm to
    // the move correct direction
    adjustment = firstDerivative / (-fabs(secondDerivative));
    // Limit the move size by .5.
    if (fabs(adjustment) < .5) {
      newEst = est - adjustment;
    } else {
      newEst = est - 0.5 * ((adjustment > 0) - (adjustment < 0));
    }
    difference = fabs(newEst - est);

    // Rcout << "  est_ability_ml_nr_iv - Stage 1.2 - FD = " <<
    //   firstDerivative << " ; SD = " << secondDerivative <<
    //     " ; Adj = " << adjustment << " ; est = " << est << " ; newEst = "
    //     << newEst << " ; diff = " << difference << std::endl;

    // If estimate is out of bounds it will bound it.
    if ((newEst <= minEst) || (newEst >= maxEst)) {

      // Rcout << "  est_ability_ml_nr_iv - Stage 1.2 - "<< difference <<
      //   std::endl;

      if ((est <= minEst) || (est >= maxEst))
      {

        // Rcout << "  est_ability_ml_nr_iv - Stage 1.3" << std::endl;

        if (est <= minEst) {
          return(minEst);
        } else {
          return(maxEst);
        }
      }
    }
    est = newEst;
  }
  return est;
}


//#############################################################################@
//########################### est_ability_4pm_nr_response_cpp ##################
//#############################################################################@
// [[Rcpp::export]]
double est_ability_4pm_nr_response_cpp(
    Rcpp::S4 resp,
    Rcpp::S4 ip,
    Rcpp::NumericVector theta_range = NumericVector::create(-5, 5),
    double criterion = 0.001,
    Rcpp::Nullable<Rcpp::NumericVector> initial_estimates = R_NilValue) {
  // Rcout << "est_ability_ml_nr - Stage 1" << std::endl;
  Rcpp::NumericVector init_est(3); // vector holding initial values
  double est; // Final estimate value that will be returned
  // Vector that holds estimates obtained using different initial values.
  bool all_equal = true;

  // resp_loglik_response_cpp functions require 'Response' objects to have
  // valid item_id slot.
  if (resp.slot("item_id") == R_NilValue)
    throw Rcpp::exception("Invalid 'resp'. 'resp' should have valid "
                            "'item_id' slot.", false);

  // Set the initial values if a NULL value provided
  if (initial_estimates.isNotNull()) {
    init_est = as<Rcpp::NumericVector>(initial_estimates);
  } else {
    init_est[0] = theta_range[0] + 2 * criterion;
    init_est[1] = 0;
    init_est[2] = theta_range[1] - 2 * criterion;
    // init_est = {theta_range[0] + 2 * criterion, 0,
    //             theta_range[1] - 2 * criterion};
  }
  int n_init_values = init_est.size();
  // Vector that holds estimates obtained using different initial values.
  Rcpp::NumericVector estimates(n_init_values);

  // Rcout << "est_ability_ml_nr - Stage 2" << std::endl;

  // Make sure that init_est has at least two values:
  if (init_est.size() < 2)
    stop("Please proivde at least two different initial values.");
  // Find ML ability estimates based on different initial values.

  // Rcpp::List ip_list = flatten_itempool_cpp(ip);
  for(int i = 0; i < n_init_values; i++) {
    estimates[i] = est_ability_4pm_nr_iv_response_cpp(
      resp, ip, theta_range, init_est[i], criterion);
    if (i > 0 && (std::fabs(estimates[i] - estimates[i-1]) >
                    n_init_values*criterion))
      all_equal = false;
    // Rcout << i << " - " << estimates[i] << std::endl;
  }
  est = estimates[0];
  // Rcout << "est_ability_ml_nr - Stage 5 - " << all_equal << " ; est = " <<
  // est << std::endl;

  // Check whether all elements of the vector are equal, if not return the
  // estimate value with the highest response log-likelihood
  if (!all_equal) {
    //Rcout << "est_ability_ml_nr - Stage 5.1 - Difference found!" << std::endl;

    // Calculate response log likelihood value at the estimate that was
    // derived from first initial estimate.
    double resp_ll = 0;
    double resp_ll_max = resp_loglik_response_cpp(estimates[0], resp, ip, 0);

    for(int i = 0; i < n_init_values; i++) {
      // Rcout << "est_ability_ml_nr - Stage 5.2" << std::endl;
      resp_ll = resp_loglik_response_cpp(estimates[i], resp, ip, 0);
      if (resp_ll > resp_ll_max) {
        est = estimates[i];
        resp_ll_max = resp_ll;
      }
    }
  }
  // Rcout << "est_ability_ml_nr - Stage 7 - " << est << std::endl;
  return(est);
}





// //##########################################################################@
// //########################### est_ability_4pm_nr_itempool_cpp ###############
// //##########################################################################@
// // [[Rcpp::export]]
// double est_ability_4pm_nr_itempool_cpp(NumericMatrix resp, Rcpp::S4 ip,
//                                         double initialEst, double criterion,
//                                         double minEst, double maxEst)
// {
//   // This function is replacement of estNR function. Does exactly the same.
//   // Estimates ability using Newton-Raphson for just one theta.
//   double est = initialEst;
//   double difference = criterion + 1;
//   NumericVector tempVector(1);
//   NumericVector tempEst(1);
//   double firstDerivative, secondDerivative, adjustment, newEst;
//
//   if ((initialEst <= minEst) || (initialEst >= maxEst)) {
//     est = 0;
//   }
//   while (difference > criterion)
//   {
//     tempEst[0] = est;
//     tempVector = resp_loglik_fd_itempool_cpp(resp, tempEst, ip);
//     firstDerivative = tempVector[0];
//
//     tempEst[0] = est;
//     tempVector = resp_loglik_sd_itempool_cpp(resp, tempEst, ip);
//     secondDerivative = tempVector[0];
//
//     // -fabs(secondDerivative) is an adjustment to make algorithm to
//     // the move correct direction
//     adjustment = firstDerivative / (-fabs(secondDerivative));
//     // Limit the move size by .5.
//     if (fabs(adjustment) < .5) {
//       newEst = est - adjustment;
//     } else {
//       newEst = est - 0.5 * ((adjustment > 0) - (adjustment < 0));
//     }
//     difference = fabs(newEst - est);
//     // If estimate is out of bounds it will bound it.
//     if ((newEst <= minEst) || (newEst >= maxEst)) {
//       if ((est <= minEst) || (est >= maxEst))
//       {
//         if (est <= minEst) {
//           return(minEst);
//         } else {
//           return(maxEst);
//         }
//       }
//     }
//     est = newEst;
//   }
//   return est;
// }


// //##########################################################################@
// //########################### est_ability_4pm_nr_matrix_cpp #################
// //##########################################################################@
// // [[Rcpp::export]]
// double est_ability_4pm_nr_matrix_cpp(
//     NumericMatrix resp, NumericMatrix ip,
//     double initialEst, double criterion, double minEst,
//     double maxEst)
// {
//   // This function is replacement of estNR function. Does exactly the same.
//   // Estimates ability using Newton-Raphson for just one theta.
//   double est = initialEst;
//   double difference = criterion + 1;
//   NumericVector tempVector(1);
//   NumericVector tempEst(1);
//   double firstDerivative, secondDerivative, adjustment, newEst;
//
//   if ((initialEst <= minEst) || (initialEst >= maxEst)) {
//     est = 0;
//   }
//   while (difference > criterion)
//   {
//     tempEst[0] = est;
//     tempVector = resp_loglik_fd_4pm_matrix_cpp(resp, tempEst, ip);
//     firstDerivative = tempVector[0];
//
//     tempEst[0] = est;
//     tempVector = resp_loglik_sd_4pm_matrix_cpp(resp, tempEst, ip);
//     secondDerivative = tempVector[0];
//
//     // -fabs(secondDerivative) is an adjustment to make algorithm to
//     // the move correct direction
//     adjustment = firstDerivative / (-fabs(secondDerivative));
//     // Limit the move size by .5.
//     if (fabs(adjustment) < .5) {
//       newEst = est - adjustment;
//     } else {
//       newEst = est - 0.5 * ((adjustment > 0) - (adjustment < 0));
//     }
//     difference = fabs(newEst - est);
//     // If estimate is out of bounds it will bound it.
//     if ((newEst <= minEst) || (newEst >= maxEst)) {
//       if ((est <= minEst) || (est >= maxEst))
//       {
//         if (est <= minEst) {
//           return(minEst);
//         } else {
//           return(maxEst);
//         }
//       }
//     }
//     est = newEst;
//   }
//   return est;
// }


//#############################################################################@
//#############################################################################@
//########################### est_ability_eap_cpp ##############################
//#############################################################################@
//#############################################################################@

//#############################################################################@
//########################### est_ability_eap_single_examinee_cpp ##############
//#############################################################################@
// [[Rcpp::export]]
Rcpp::List est_ability_eap_single_examinee_cpp(Rcpp::NumericVector resp,
                                               Rcpp::S4 ip,
                                               Rcpp::NumericVector theta_range,
                                               int no_of_quadrature,
                                               Rcpp::String prior_dist,
                                               Rcpp::NumericVector prior_par)
{
  // Rcout << "  (est_ability_eap_single_examinee_cpp) -- 1 -- Beginning" <<
  // std::endl;

  Rcpp::List item_list = as<List>(ip.slot("item_list"));
  unsigned int noItem = item_list.size();
  double L; // This will hold the likelihood of a response string
  Rcpp::NumericVector x(no_of_quadrature);
  Rcpp::NumericVector fx_numerator(no_of_quadrature);
  Rcpp::NumericVector fx_denominator(no_of_quadrature);
  Rcpp::NumericVector fx_std_error(no_of_quadrature);
  double prior; // The value of the prior distribution
  // this will hold numerator of the posterior mean of theta
  double posterior_numerator;
  double posterior_denominator; // this will hold denominator of the posterior
  double est;
  double se;
  S4 temp_item("Item");
  Rcpp::List output;

  // Calculate the theta point vector
  for (int i = 0; i < no_of_quadrature; i++)
    x[i] = theta_range[0] + i * (theta_range[1] -
      theta_range[0]) / (no_of_quadrature - 1);
  for (int i = 0; i < no_of_quadrature; i++)
  {
    // Calculate the likelihood
    L = 1;
    for (unsigned int j = 0; j < noItem; j++)
    {
      // Only use non-missing, i.e. non-NA items.
      if (!NumericVector::is_na(resp[j])) {
        L = L * resp_lik_bare_item_cpp(resp[j], x[i],
                                       as<Rcpp::S4>(item_list[j]));
      }
    }
    // Rcout << "      (est_ability_eap_single_examinee_cpp) -- 3.1.3" <<
    // std::endl;
    // Calculate the prior distribution value
    if (prior_dist == "norm") {
      // mean = prior_par[0], sd = prior_par[1]
      prior = R::dnorm(x[i], prior_par[0], prior_par[1], false);
    } else if (prior_dist == "unif") {
      // min = prior_par[0], max = prior_par[1], log = false
      prior = R::dunif(x[i], prior_par[0], prior_par[1], false);
    } else if (prior_dist == "lnorm") {
      // meanlog = prior_par[0], sdlog = prior_par[1], log = false
      prior = R::dlnorm(x[i], prior_par[0], prior_par[1], false);
    } else if (prior_dist == "gamma") {
      // shape = prior_par[0], scale = prior_par[1], log = false
      prior = R::dgamma(x[i], prior_par[0], prior_par[1], false);
    } else if (prior_dist == "t") {
      // df = prior_par[0], log = false
      prior = R::dt(x[i], prior_par[0], false);
    } else if (prior_dist == "cauchy") {
      // location = prior_par[0], scale = prior_par[1], log = false
      prior = R::dcauchy(x[i], prior_par[0], prior_par[1], false);
    } else {
      // Todo: Double check this default value
      prior = 1;
    }
    fx_numerator[i] = L * prior * x[i];
    fx_denominator[i] = L * prior;
  }

  posterior_numerator = integrate(x, fx_numerator);
  posterior_denominator = integrate(x, fx_denominator);
  est = posterior_numerator/posterior_denominator;
  for (int i = 0; i < no_of_quadrature; i++)
    fx_std_error[i] = (x[i] - est) * (x[i] - est) * fx_denominator[i];
  se = sqrt(integrate(x, fx_std_error) / posterior_denominator);
  // Rprintf("r: %d; est: %.3f\n", r, est[r]);
  // Estimate standard error
  output["est"] = est;
  output["se"] = se;
  return output;
}


//#############################################################################@
//########################### est_ability_eap_cpp ##############################
//#############################################################################@
// Estimate ability for a group of subjects, in 'resp' rows represents
// examinees, and columns represent items.
// [[Rcpp::export]]
Rcpp::List est_ability_eap_cpp(Rcpp::NumericMatrix resp,
                               Rcpp::S4 ip,
                               Rcpp::NumericVector theta_range,
                               int no_of_quadrature,
                               Rcpp::String prior_dist,
                               Rcpp::NumericVector prior_par)
{
  int no_of_examinees = resp.nrow();
  Rcpp::NumericVector est(no_of_examinees);
  Rcpp::NumericVector se(no_of_examinees);
  Rcpp::List output, temp_output;
  Rcpp::NumericVector temp_resp;

  for (int r = 0; r < no_of_examinees; r++) {
    temp_resp = resp(r, _);
    temp_output = est_ability_eap_single_examinee_cpp(
      temp_resp, ip, theta_range, no_of_quadrature, prior_dist, prior_par);
    est[r] = temp_output["est"];
    se[r] = temp_output["se"];
  }
  // Estimate standard error
  output["est"] = est;
  output["se"] = se;
  return output;
}



//#############################################################################@
//########################### est_ability_eap_response_cpp #####################
//#############################################################################@
//' @param resp A Response object.
//' @param ip An Itempool object.
//' @noRd

Rcpp::List est_ability_eap_response_cpp(
    Rcpp::S4 resp,
    Rcpp::S4 ip,
    Rcpp::NumericVector theta_range = NumericVector::create(-5, 5),
    int no_of_quadrature = 61,
    std::string prior_dist = "norm",
    Rcpp::NumericVector prior_par = NumericVector::create(0, 1))
{
  double L; // This will hold the likelihood of a response string
  Rcpp::NumericVector x(no_of_quadrature);
  Rcpp::NumericVector fx_numerator(no_of_quadrature);
  Rcpp::NumericVector fx_denominator(no_of_quadrature);
  Rcpp::NumericVector fx_std_error(no_of_quadrature);
  double prior; // The value of the prior distribution
  // this will hold numerator of the posterior mean of theta
  double posterior_numerator;
  double posterior_denominator; // this will hold denominator of the posterior
  double est;
  double se;
  Rcpp::S4 temp_item("Item");
  Rcpp::List output;

  // Calculate the theta point vector
  for (int i = 0; i < no_of_quadrature; i++)
    x[i] = theta_range[0] + i * (theta_range[1] -
      theta_range[0]) / (no_of_quadrature - 1);

  for (int i = 0; i < no_of_quadrature; i++) {
    // Calculate the likelihood
    L = resp_lik_response_cpp(x[i], resp, ip);
    // Calculate the prior distribution value
    if (prior_dist == "norm") {
      // mean = prior_par[0], sd = prior_par[1]
      prior = R::dnorm(x[i], prior_par[0], prior_par[1], false);
    } else if (prior_dist == "unif") {
      // min = prior_par[0], max = prior_par[1], log = false
      prior = R::dunif(x[i], prior_par[0], prior_par[1], false);
    } else if (prior_dist == "lnorm") {
      // meanlog = prior_par[0], sdlog = prior_par[1], log = false
      prior = R::dlnorm(x[i], prior_par[0], prior_par[1], false);
    } else if (prior_dist == "gamma") {
      // shape = prior_par[0], scale = prior_par[1], log = false
      prior = R::dgamma(x[i], prior_par[0], prior_par[1], false);
    } else if (prior_dist == "t") {
      // df = prior_par[0], log = false
      prior = R::dt(x[i], prior_par[0], false);
    } else if (prior_dist == "cauchy") {
      // location = prior_par[0], scale = prior_par[1], log = false
      prior = R::dcauchy(x[i], prior_par[0], prior_par[1], false);
    } else {
      // Todo: Double check this default value
      prior = 1;
    }
    fx_numerator[i] = L * prior * x[i];
    fx_denominator[i] = L * prior;
  }

  posterior_numerator = integrate(x, fx_numerator);
  posterior_denominator = integrate(x, fx_denominator);
  est = posterior_numerator/posterior_denominator;
  // Estimate standard error
  for (int i = 0; i < no_of_quadrature; i++)
    fx_std_error[i] = (x[i] - est) * (x[i] - est) * fx_denominator[i];
  se = sqrt(integrate(x, fx_std_error) / posterior_denominator);
  output["est"] = est;
  output["se"] = se;
  return output;
}


//#############################################################################@
//########################### est_ability_eap_response_set_cpp #################
//#############################################################################@
//' @param resp_set A Response_set object.
//' @param ip An Itempool object.
//' @noRd

// [[Rcpp::export]]
Rcpp::List est_ability_eap_response_set_cpp(
    Rcpp::S4 resp_set,
    Rcpp::S4 ip,
    Rcpp::NumericVector theta_range = NumericVector::create(-5, 5),
    int no_of_quadrature = 61,
    std::string prior_dist = "norm",
    Rcpp::NumericVector prior_par = NumericVector::create(0, 1)) {

  // Check whether resp_set is Response_set object. If it is a Response object
  // convert it to a Response_set object with one Response.
  if (resp_set.inherits("Response")) {
    Rcpp::S4 response_set("Response_set");
    Rcpp::List response_list(1);
    response_list[0] = resp_set;
    if (resp_set.slot("examinee_id") == R_NilValue) {
      response_list.attr("names") = "Ex-1";
    } else {
      response_list.attr("names") = as<Rcpp::StringVector>(
        resp_set.slot("examinee_id"));
    }
    response_set.slot("response_list") = response_list;
    resp_set = response_set;
  } else if (!resp_set.inherits("Response_set")) {
    stop("Invalid 'resp_set' argument. 'resp_set' should be a Response_set "
           "object.");
  }
  // Make sure resp_set and ip are valid and compatible
  check_validity_response_set_cpp(resp_set, ip);
  // Rcpp::List ip_list = flatten_itempool_cpp(ip);

  Rcpp::List resp_list = as<Rcpp::List>(resp_set.slot("response_list"));
  int num_of_resp = resp_list.size();

  Rcpp::NumericVector est(num_of_resp);
  Rcpp::NumericVector se(num_of_resp);
  Rcpp::S4 temp_resp;
  Rcpp::List output, temp_output;
  for (int i = 0; i < num_of_resp; i++) {
    temp_resp = as<Rcpp::S4>(resp_list[i]);
    temp_output = est_ability_eap_response_cpp(
      temp_resp, ip, theta_range, no_of_quadrature, prior_dist, prior_par);
    est[i] = temp_output["est"];
    se[i] = temp_output["se"];
  }
  // Estimate standard error
  output["est"] = est;
  output["se"] = se;
  return output;
}

//#############################################################################@
//#############################################################################@
//########################### est_ability_owen #################################
//#############################################################################@
//#############################################################################@

//#############################################################################@
//########################### est_ability_owen_item_cpp ########################
//#############################################################################@
// [[Rcpp::export]]
Rcpp::List est_ability_owen_item_cpp(Rcpp::S4 item, int resp,
                                     double m0, double v0) {

  // Item difficulty
  double b = as<double>(item.slot("b"));
  // Item discrimination
  double D = item.hasSlot("D") ? as<double>(item.slot("D")) : 1;
  double a = item.hasSlot("a") ? as<double>(item.slot("a")) : 1;
  a = a/D;
  double c = item.hasSlot("c") ? as<double>(item.slot("c")) : 0;
  double d = item.hasSlot("d") ? as<double>(item.slot("d")) : 1;
  // DD is the D in original owen's article. Here I rename it in case of
  // a confusion with scaling parameter D, 1.7.
  double DD = (b - m0) / sqrt((1/(a*a)) + v0);
  double dnormDD, pnormDD;
  NumericVector temp_nv;
  NumericVector DD_vector = NumericVector::create(DD);
  temp_nv = pnorm(DD_vector);
  pnormDD = temp_nv[0];
  temp_nv = dnorm(DD_vector);
  dnormDD = temp_nv[0];
  DD_vector[0] = -DD;
  temp_nv = pnorm(DD_vector);
  double pnormDD_minus = temp_nv[0];
  double A = c + (d - c) * pnormDD_minus;
  double m1 = m0 - v0 * (1/sqrt((1/(a*a)) + v0)) * (dnormDD / pnormDD) * (
    1 - resp / A);
  // double m1 = m0 - v0 * pow(pow(a, -2) + v0, -.5) * (dnormDD / pnormDD) * (
  //   1 - resp / A);
  double v1 = v0 - v0*v0 * (1/((1/(a*a)) + v0)) * (
    dnormDD / pnormDD) * (1 - resp / A) * ((dnormDD / pnormDD) * (
        1 - resp / A) + DD);
  // double v1 = v0 - pow(v0, 2) * pow(pow(a, -2) + v0, -1) * (
  //   dnormDD / pnormDD) * (1 - resp / A) * ((dnormDD / pnormDD) * (
  //       1 - resp / A) + DD);
  List output;
  output["m1"] = m1;
  output["v1"] = v1;
  return output;
}


//#############################################################################@
//########################### est_ability_owen_cpp #############################
//#############################################################################@
// [[Rcpp::export]]
Rcpp::List est_ability_owen_cpp(Rcpp::S4 ip, Rcpp::NumericVector resp,
                                double m0, double v0) {
  // Rcout << "Stage (EstAbilityOwen) 1" << std::endl;
  List output;
  // if (Rf_inherits(ip, "Item")) {
  if (Rf_inherits(ip, "Rasch") || Rf_inherits(ip, "1PL") ||
      Rf_inherits(ip, "2PL") || Rf_inherits(ip, "3PL") ||
      Rf_inherits(ip, "4PL")) {
    resp = resp[0];
    output = est_ability_owen_item_cpp(ip, as<int>(resp), m0, v0);
  } else if (Rf_inherits(ip, "Itempool")) {
    // Rcout << "Stage (EstAbilityOwen) 2" << std::endl;
    int no_items = resp.size();
    // Rcout << "Stage (EstAbilityOwen) 3" << std::endl;
    List item_list = ip.slot("item_list");
    output = List::create(Named("m1") = m0, Named("v1") = v0);
    // Rcout << "Stage (EstAbilityOwen) 4" << std::endl;
    for (int i = 0; i < no_items; i++) {
      if (!NumericVector::is_na(resp[i]))
        output = est_ability_owen_item_cpp(as<Rcpp::S4>(item_list(i)), resp[i],
                                           output["m1"], output["v1"]);
    }
    // Rcout << "Stage (EstAbilityOwen) 5" << std::endl;
  } else
    stop("This function can only process an 'Item' and 'Itempool' objects.");
  // The final output is given as the square root of the posterior variance v1.
  return List::create(Named("est") = output["m1"],
                      Named("se") = sqrt(as<double>(output["v1"])));
}



//#############################################################################@
//#############################################################################@
//########################### est_ability_optim ################################
//#############################################################################@
//#############################################################################@

//#############################################################################@
//########################### est_ability_optim_response_cpp ###################
//#############################################################################@
//' @param resp A Response object.
//' @param ip_list A list of Item/Testlet objects.
//' @noRd

// [[Rcpp::export]]
double est_ability_optim_response_cpp(
    Rcpp::S4 resp,
    Rcpp::List ip_list,
    Rcpp::NumericVector theta_range = NumericVector::create(-5, 5),
    double tol = 0.0000001) {

  // Extract R's optim function
  Rcpp::Environment stats("package:stats");
  Rcpp::Function optim = stats["optim"];

  // Call the optim function from R in C++
  Rcpp::List opt_results = optim(
    Rcpp::_["par"]    = Rcpp::NumericVector::create(0),
    // Make sure this function is not exported!
    Rcpp::_["fn"]     = Rcpp::InternalFunction(&resp_loglik_response_cpp),
    Rcpp::_["method"] = "Brent",
    Rcpp::_["lower"] = Rcpp::NumericVector::create(-5),
    Rcpp::_["upper"] = Rcpp::NumericVector::create(5),
    // Pass in the other parameters as everything
    // is scoped environmentally
    Rcpp::_["resp"] = resp,
    Rcpp::_["ip_list"] = ip_list,
    Rcpp::_["derivative"] = 0
    );

  return Rcpp::as<double>(opt_results[0]);
}


