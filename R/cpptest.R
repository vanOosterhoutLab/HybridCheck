src <- '
Rcpp::CharacterVector mystring(instring);
Rcpp::NumericVector out(1); 
out[0] = mystring.length();
return out;
'

fun <- cxxfunction(signature(instring="character"), src, plugin="Rcpp")