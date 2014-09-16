#include <Rcpp.h>
using namespace Rcpp;

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

DataFrame ScanLoop_cpp(DataFrame df, CharacterMatrix seq){
  std::vector<int> frompos = Rcpp::as< std::vector<int> >(df["WindowStart"]);
  std::vector<int> topos = Rcpp::as< std::vector<int> >(df["WindowEnd"]);
  std::vector<double> AB = Rcpp::as< std::vector<double> >(df["AB"]);
  std::vector<double> AC = Rcpp::as< std::vector<double> >(df["AC"]);
  std::vector<double> BC = Rcpp::as< std::vector<double> >(df["BC"]);
  
  for (int i = 0; i < AB.size(); i++){
    CharacterMatrix dnaStretch = seq(_, Range(frompos[i], topos[i]));
    
  }
  
}
