Foo <- setRefClass("Foo",
                   fields = list(
                     file = "character",
                     x = function(value) {
                       if(missing(value)){
                         read.table( file )
                       } else {
                         write.table(value, file = file)
                       }
                     }
                     ),
                   methods = list(
                     initialize =
                       function( path ) {
                         file <<- path
                         x <<- data.frame(A = NA, B = NA, C = NA)
                       })
                   )

testfoo <- Foo$new(file="~/Desktop/testfoo")

FooFunc <- function( x ){
  x$x[ 2, ] <- c(7,8,9)
  x$x[3,1] <- 2
}
