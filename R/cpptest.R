src <- '
//Convert the inputted character matrix of DNA sequences an Rcpp class.
Rcpp::CharacterMatrix mymatrix(inmatrix);

//Get the number of columns and rows in the matrix
int ncolumns = mymatrix.ncol();
int numrows = mymatrix.nrow();

//Get the dimension names
Rcpp::List dimnames = mymatrix.attr("dimnames");

Rcpp::CharacterMatrix vec1 = mymatrix(Range(1,numrows),_);
Rcpp::CharacterMatrix vec2 = mymatrix(Range(0,numrows-1),_); 
'
myCut <- function(FullMatrix){
  return( FullMatrix[, colSums( FullMatrix[-1,] != FullMatrix[-nrow( FullMatrix ), ] ) > 0] )
}

seqCleaner <- cxxfunction(signature(inmatrix="character"), src, plugin="Rcpp")

Rcpp::NumericMatrix Am(A);
+     int nrows = Am.nrow();
+     int ncolumns = Am.ncol();
+     for (int i = 0; i < ncolumns; i++) {
  +         for (int j = 1; j < nrows; j++) {
    +             Am(j,i) = Am(j,i) + Am(j-1,i);
    +         }
  +     }
+     return Am;

List dimnames = x.attr( "dimnames" ) ;
return dimnames[0] ;



src <- '
#include<vector>
CharacterMatrix inDNA(input);
std::vector<int> informativeSites; 
  for(int i = 0; i < inDNA.ncol(); i++)
  {
    CharacterMatrix bpsite = inDNA(_,i);
    if(all(bpsite == bpsite[1]))
    {
      informativeSites.push_back(i);
    }
  }
CharacterMatrix cutDNA = completeDNA(,informativeSites);
return cutDNA;
'

require( Rcpp )
cppFunction('
  LogicalVector unq_mat( CharacterMatrix x ){

  int nc = x.ncol() ;
  LogicalVector out(nc);

  for( int i=0; i < nc; i++ ) {
    out[i] = unique( x(_,i) ).size() != 1 ;
    }
  return out;
}'
)

cppFunction('
  LogicalVector findN( CharacterMatrix x ){

  int nc = x.ncol() ;
  LogicalVector out(nc);
  CharacterVector ns(x.nrow());
  std::fill(ns.begin(), ns.end(), "n");
  

  for( int i=0; i < nc; i++ ) {
    Rcpp::CharacterVector vec = x(_,i); 
    out[i] = is_true(any(vec == ns)) ;
    }
  return out;
}'
)


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
  return(input[,-unq_mat(input)])
}

mat <- matrix( as.character(c(rep(1,5),sample(3,15,repl=TRUE),rep(5,5))),5)