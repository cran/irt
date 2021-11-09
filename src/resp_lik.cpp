#include <Rcpp.h>
#include "misc.h"
#include "itempool_class_methods.h"
#include "response_set_methods.h"
#include "prob.h"
using namespace Rcpp;

//#############################################################################@
//#############################################################################@
//########################### resp_lik #########################################
//#############################################################################@
//#############################################################################@


//#############################################################################@
//########################### resp_lik_bare_item_cpp ###########################
//#############################################################################@
// This function calculates the response likelihood of an item for a
// given response.

// [[Rcpp::export]]
double resp_lik_bare_item_cpp(double resp, double theta, Rcpp::S4 item) {
  // Deal with missing responses, return NA directly
  if (NumericVector::is_na(resp))
    return NA_REAL;

  // Get the Psychometric Model name
  std::string model = as<std::string>(item.attr("class"));

  if (model == "GPCM" || model == "GPCM2" || model == "PCM" || model == "GRM") {
    Rcpp::NumericVector P;
    if (model == "GPCM" || model == "PCM" || model == "GPCM2") {
      P = prob_gpcm_bare_cpp(theta, item, 0, resp);
      return P[0];
    } else if (model == "GRM") {
      P = prob_grm_bare_cpp(theta, item);
      return P[resp];
    }
    // The following line assumes that the resp goes from 0 to maximum number
    // of categories.
    return P[resp];
  } else if (check_item_model(item, true, true)) {
    return prob_4pm_bare_cpp(theta, item, 0, resp);
    // // The following line is important (instead of second line) because
    // // it accounts for resp values that are not 0 or 1.
    // return pow(P, resp) * pow(1.0-P, 1.0-resp);
    // return resp * P + (1 - resp) * (1 - P);
  }
  return NA_REAL;
}


//##############################################################################
//########################### resp_lik_item_cpp ################################
//##############################################################################
// [[Rcpp::export]]
Rcpp::NumericVector resp_lik_item_cpp(Rcpp::NumericVector resp,
                                      Rcpp::NumericVector theta, Rcpp::S4 item)
{
  // Calculate response likelihood for one item and multiple theta's (and
  // responses)
  unsigned int num_of_theta = theta.size();
  Rcpp::NumericVector output(num_of_theta);
  for(unsigned int i = 0; i < num_of_theta; i++)
    output[i] = resp_lik_bare_item_cpp(resp[i], theta[i], item);
  return output;
}


//##############################################################################
//########################### resp_lik_bare_testlet_cpp ########################
//##############################################################################
// [[Rcpp::export]]
double resp_lik_bare_testlet_cpp(Rcpp::NumericVector resp, double theta,
                                 Rcpp::S4 testlet)
{
  // Calculate response log-likelihood for a testlet and single examinee
  Rcpp::List item_list = as<List>(testlet.slot("item_list"));
  unsigned int num_of_items = item_list.size();
  double output = 1;
  Rcpp::S4 item; // This will be item
  for(unsigned int i = 0; i < num_of_items; i++) {
    item = as<Rcpp::S4>(item_list(i));
    if (!NumericVector::is_na(resp[i]))
        output = output * resp_lik_bare_item_cpp(
          resp(i), theta, as<Rcpp::S4>(item_list(i)));
  }
  return output;
}


//##############################################################################
//########################### resp_lik_testlet_cpp #############################
//##############################################################################
// [[Rcpp::export]]
Rcpp::NumericVector resp_lik_testlet_cpp(Rcpp::NumericMatrix resp,
                                         Rcpp::NumericVector theta,
                                         Rcpp::S4 testlet)
{
  // Calculate response log-likelihood for an Itempool and multiple
  // theta's (and response strings)
  unsigned int num_of_theta = theta.size();
  NumericVector output(num_of_theta);
  for(unsigned int i = 0; i < num_of_theta; i++) {
    // Get the row belong to the examinee. It is assumed that each row represents
    // an examinee.
    NumericVector resp_vector = resp(i, _);
    output[i] = resp_lik_bare_testlet_cpp(resp_vector, theta[i], testlet);
  }
  return output;
}


//#############################################################################@
//########################### resp_lik_bare_itempool_cpp ######################
//#############################################################################@

// [[Rcpp::export]]
double resp_lik_bare_itempool_cpp(Rcpp::NumericVector resp, double theta,
                                   Rcpp::S4 ip) {
  // This function calculates the response likelihood of an item pool for a
  // given response string and one theta.

  // Assuming that theta.size() ==  resp.size(), though it may not be the case
  // always. There might be an instance where more responses can happen if
  // item pool has testlets.
  int no_items = resp.size();
  double result = 1;
  S4 item;
  // Indicator variable for whether all responses are missing (true) or at
  // least there is one non-missing response (false).
  bool resp_all_na = true;
  List item_list = ip.slot("item_list");
  for (int i = 0; i < no_items; i++) {
    // iterate over non-missing responses
    if (!R_IsNA(resp[i])) {
      resp_all_na = false; // one non-na observed
      item = as<Rcpp::S4>(item_list[i]);
      result = result * resp_lik_bare_item_cpp(resp[i], theta, item);
    }
  }
  if (resp_all_na) result = NA_REAL;   // should it return 0 or NA?
  return result;
}


//##############################################################################
//########################### resp_lik_itempool_cpp ###########################
//##############################################################################
// [[Rcpp::export]]
Rcpp::NumericVector resp_lik_itempool_cpp(Rcpp::NumericMatrix resp,
                                          Rcpp::NumericVector theta,
                                          Rcpp::S4 ip)
{
  // Calculate response log-likelihood for an Itempool and multiple
  // theta's (and response strings)
  unsigned int num_of_theta = theta.size();
  NumericVector output(num_of_theta);
  for(unsigned int i = 0; i < num_of_theta; i++) {
    // Get the row belong to the examinee. It is assumed that each row represents
    // an examinee.
    NumericVector resp_vector = resp(i, _);
    output[i] = resp_lik_bare_itempool_cpp(resp_vector, theta[i], ip);
  }
  return output;
}




//#############################################################################@
//########################### resp_lik_response_cpp ############################
//#############################################################################@
// @param resp A Response object.


// [[Rcpp::export]]
double resp_lik_response_cpp(double theta, Rcpp::S4 resp, Rcpp::S4 ip)
{
  Rcpp::NumericVector scores = as<Rcpp::NumericVector>(resp.slot("score"));
  Rcpp::StringVector item_ids = as<Rcpp::StringVector>(resp.slot("item_id"));

  double output = 1;
  Rcpp::S4 item; // This will be an item or testlet
  std::string item_id;
  Rcpp::List ip_list = as<Rcpp::List>(ip.slot("item_list"));
  int num_of_items = scores.size();

  // check if testlet_id is NULL, TYPEOF(R_NilValue)
  if (TYPEOF(resp.slot("testlet_id")) == 0) {
    for (int i = 0; i < num_of_items; i++) {
      item_id = item_ids[i];
      item = as<Rcpp::S4>(ip_list[item_id]);
      output = output * resp_lik_bare_item_cpp(scores[i], theta, item);
    }
  } else {
    int item_no = 0;
    std::string testlet_id;
    // The vector holding testlet_ids column of resp.
    Rcpp::StringVector testlet_ids = as<Rcpp::StringVector>(
      resp.slot("testlet_id"));
    Rcpp::StringVector calculated_teslelet_ids;
    Rcpp::NumericVector selected_testlet_item_scores;
    Rcpp::S4 temp_s4;

    // holds whether a testlet administered before or not
    bool administered = false;
    while (item_no < num_of_items) { // Standalone item
      item_id = item_ids[item_no];
      if (Rcpp::StringVector::is_na(testlet_ids[item_no])) {
        item = as<Rcpp::S4>(ip_list[item_id]);
        output = output * resp_lik_bare_item_cpp(scores[item_no], theta, item);
      } else { // Testlet item
        testlet_id = testlet_ids[item_no];
        // Check if the likelihood of this testlet has already been calculated
        administered = false;
        for (int j = 0; j < calculated_teslelet_ids.size(); j++) {
          if (testlet_id == as<std::string>(calculated_teslelet_ids[j]))
            administered = true;
        }
        // Find all of the testlet items.
        if (!administered) { // Calculate the likelihood of all items in the
                             // testlet
          // Get testlet object:
          item = as<Rcpp::S4>(ip_list[testlet_id]);
          // item pool of the testlet
          temp_s4 = as<Rcpp::S4>(item.slot("item_list"));
          Rcpp::NumericVector selected_testlet_item_scores(
              get_itempool_size(temp_s4)[0], Rcpp::NumericVector::get_na());

          selected_testlet_item_scores.attr("names") = get_ids_itempool_cpp(
            temp_s4);
          // Get the scores the testlet items in the same order that appears
          // in item pool. If there is an unadministered item, let it be NA
          // testlet_scores
          for (int j = item_no; j < num_of_items; j++) {
            if (testlet_ids[j] ==  testlet_id) { // if item belongs to this
                                                     // testlet..
              item_id = item_ids[j];
              selected_testlet_item_scores[item_id] = scores[j];
            }
          }
          // Calculate likelihood of the testlet
          output = output * resp_lik_bare_testlet_cpp(
            selected_testlet_item_scores, theta, item);
          calculated_teslelet_ids.push_back(testlet_id);
        } // else -> likelihood of testlet item has already been calculated
      }
      item_no++;
    }
  }
  return output;
}



// double resp_lik_response_cpp(Rcpp::S4 resp, double theta, Rcpp::List ip_list) {
//   Rcpp::NumericVector scores = as<Rcpp::NumericVector>(resp.slot("score"));
//   Rcpp::StringVector item_ids = as<Rcpp::StringVector>(resp.slot("item_id"));
//
//   int num_of_items = scores.size();
//
//   double output = 1;
//   Rcpp::S4 item; // This will be item
//   std::string item_id;
//   for (int i = 0; i < num_of_items; i++) {
//     item_id = item_ids[i];
//     item = as<Rcpp::S4>(ip_list[item_id]);
//     output = output * resp_lik_bare_item_cpp(scores[i], theta, item);
//   }
//   return output;
// }


//#############################################################################@
//########################### resp_lik_response_set_cpp ########################
//#############################################################################@
// @param resp_set A Response_set object.

// [[Rcpp::export]]
Rcpp::NumericVector resp_lik_response_set_cpp(
    Rcpp::S4 resp_set,
    Rcpp::NumericVector theta,
    Rcpp::S4 ip) {

  // Make sure resp_set and ip are valid and compatible
  check_validity_response_set_cpp(resp_set, ip);
  // Rcpp::List ip_list = flatten_itempool_cpp(ip);

  Rcpp::List resp_list = as<Rcpp::List>(resp_set.slot("response_list"));
  int num_of_resp = resp_list.size();

  if (theta.size() != num_of_resp)
    stop("Incompatible 'theta' and 'resp_set'. Their length should be equal.");

  Rcpp::NumericVector output(num_of_resp);
  Rcpp::S4 temp_resp;
  for (int i = 0; i < num_of_resp; i++) {
    temp_resp = as<Rcpp::S4>(resp_list[i]);
    output[i] = resp_lik_response_cpp(theta[i], temp_resp, ip);
  }
  return output;
}






