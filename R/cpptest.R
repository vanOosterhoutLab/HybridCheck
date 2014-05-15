require(Rcpp)

cppFunction('
  LogicalVector unq_mat( CharacterMatrix x ){
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
}'
)

cppalt <- function(input){
  return(input[,unq_mat(input)])
}

mat <- matrix(as.character(c(rep(1,5),sample(3,15,repl=TRUE),rep(5,5))),5)