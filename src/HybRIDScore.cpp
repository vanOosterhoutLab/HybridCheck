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
  std::vector<int> frompos = as< std::vector<int> >(df["WindowStart"]);
  std::vector<int> topos = as< std::vector<int> >(df["WindowEnd"]);
  NumericVector AB = df["AB"];
  NumericVector AC = df["AC"];
  NumericVector BC = df["BC"];
  
  for (int i = 0; i < AB.size(); i++){
    CharacterMatrix dnaStretch = seq(_, Range(frompos[i], topos[i]));
    AB[i] = (sum(dnaStretch(1, _) == dnaStretch(2, _)) / dnaStretch.ncol()) * 100;
    AC[i] = (sum(dnaStretch(1, _) == dnaStretch(3, _)) / dnaStretch.ncol()) * 100;
    BC[i] = (sum(dnaStretch(2, _) == dnaStretch(3, _)) / dnaStretch.ncol()) * 100;
  }
  
  
}
