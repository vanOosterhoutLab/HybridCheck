#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
LogicalVector sequenceChecker_cpp( CharacterMatrix x ){
  int nc = x.ncol() ;
  CharacterVector ns(x.nrow());
  CharacterVector Ns(x.nrow());
  std::fill(ns.begin(), ns.end(), "n");
  std::fill(Ns.begin(), Ns.end(), "N");
  LogicalVector out(nc);
  for( int i=0; i < nc; i++ ) {
    out[i] = is_false(any(x(_,i) == Ns)) && is_false(any(x(_,i) == ns)) && unique( x(_,i) ).size() != 1 ;
    }
  return out;
}
