Foo <- setRefClass("Foo",
                   fields = list(
                     file = "character",
                     x = function(value) {
                       if(missing(value)){
                         read.table( file )
                       } else {
                         write.table(value, file = file)
                       }
                     },
                     y = function(value) {
                       if(missing(value)){
                         read.table( file )
                       } else {
                         write.table(value, file=file)
                       }
                     }
                     )
                   )



