#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
SEXP helloEcpp(SEXP greetingPointer) {
  Rcpp::StringVector greeting(greetingPointer);
  Rcpp::NumericVector result(greeting.size());
  for (int i=0; i<greeting.size(); i++) {
    result[i] = strlen(greeting[i]);
  }
  return(result);
}
