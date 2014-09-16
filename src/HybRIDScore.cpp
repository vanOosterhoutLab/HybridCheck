#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
LogicalVector sequenceChecker_cpp(CharacterMatrix x){
  int nc = x.ncol() ;
  CharacterVector ns(x.nrow());
  CharacterVector Ns(x.nrow());
  std::fill(ns.begin(), ns.end(), "n");
  std::fill(Ns.begin(), Ns.end(), "N");
  LogicalVector out(nc);
  for(int i=0; i < nc; i++){
    out[i] = is_false(any(x(_,i) == Ns)) && is_false(any(x(_,i) == ns)) && unique(x(_,i)).size() != 1;
  }
  return out;
}

DataFrame ScanLoop_cpp(DataFrame df){
  std::vector<double> AB = Rcpp::as< std::vector<double> >(df["AB"]);
  std::vector<double> AC = Rcpp::as< std::vector<double> >(df["AC"]);
  std::vector<double> BC = Rcpp::as< std::vector<double> >(df["BC"]);
  
  for (int i = 0; i < AB.size() ; i++){
    v_prob[i] = vaccinate_cxx(age_v[i],female_v[i],ily_v[i]); 
  }
  
}
