# Iterators over and across sequence data.

setOldClass(c("abstractiter", "iter"))

Iterator <- setRefClass(
  "Iterator",
  fields = list(checkFun = "function"),
  methods = list(
    initialize = function(fun = function()
      TRUE) {
      checkFun <<- fun
    },
    nextElem = function()
      stop("Not implemented")
  )
)

IdxIterator <- setRefClass(
  "IdxIterator",
  fields = list(current = "integer",
                limit = "integer"),
  contains = "Iterator",
  methods = list(
    initialize = function(limit = 1L, fun = function(...)
      TRUE) {
      callSuper(fun)
      current <<- 1L
      limit <<- limit
    }
  )
)

PairIdx <- setRefClass(
  "PairIdx",
  contains = "IdxIterator",
  methods = list(
    initialize = function(limit = 1L, fun = function(...)
      TRUE) {
      callSuper(limit, fun)
      current <<- c(1L, 1L)
    },
    nextElem = function() {
      repeat {
        if (!hasNext()) {
          stop("StopIteration")
        }
        if (current[2] == limit) {
          current[1] <<- current[1] + 1L
          current[2] <<- current[1]
        } else {
          current[2] <<- current[2] + 1L
        }
        if (checkFun(current)) {
          return(current)
        }
      }
    },
    hasNext = function() {
      return((current[1] <= limit) && (current[2] <= limit))
    }
  )
)

OneVAll <- setRefClass("OneVAll",
                       fields = list(one = "integer"),
                       contains = "PairIdx",
                       methods = list(
                         initialize = function(limit = 1L, one = 1L, fun = function(...) TRUE){
                           callSuper(limit, fun)
                           one <<- one
                         },
                         nextElem = function() {
                           repeat {
                             if (!hasNext()) {
                               stop("StopIteration")
                             }
                             current[2] <<- current[2] + 1L
                             if (checkFun(current)) {
                               return(current)
                             }
                           }
                         },
                         hasNext = function(){
                           return(current[2] != limit)
                         }
                       )
)
                       
                       

setRefClass("SequenceIterator",
            contains = c("VIRTUAL"))

PairsMAlign <- setRefClass(
  "PairsMAlign",
  fields = list(obj = "MultipleAlignment",
                itr = "PairIdx",
                checkFun = "function"),
  contains = "SequenceIterator",
  methods = list(
    initialize = function(obj, fun = function(x)
      TRUE) {
      obj <<- obj
      itr <<- PairIdx$new(limit = nrow(obj))
    },
    nextElem = function() {
      repeat {
        idx <- itr$nextElem()
        pair <- maskSequences(obj,
                              idx,
                              invert = TRUE,
                              append = "replace")
        if(checkFun(pair)){
          return(pair)
        }
      }
    }
  )
)

RefMAlign <- setRefClass(
  "RefMAlign",
  fields = list(
    obj = "MultipleAlignment",
    itr = "OneVAll",
    ref = "integer"),
  contains = "SequenceIterator",
  methods = list(
    initialize = function(obj, ref = NULL, fun = function(x)
      TRUE) {
      if(is.null(ref)){
        ref <<- 1L
      } else {
        if(class(ref) == "character"){
          ref <<- which(rownames(obj)) == ref
        } else {
          if(class(ref) == "integer" || class(ref) == "numeric"){
            ref <<- ref
          } else {
            stop("Invalid input for parameter 'ref'.")
          }
        }
      }
      itr <<- OneVAll$new(limit = nrow(obj), one = ref)
    },
    nextElem = function() {
      repeat {
        idx <- itr$nextElem()
        pair <- maskSequences(obj,
                              idx,
                              invert = TRUE,
                              append = "replace")
        if(checkFun(pair)){
          return(pair)
        }
      }
    }
  )
)

setGeneric("Pairs", function(object){
  standardGeneric("Pairs")
})

setMethod("Pairs",
          representation("MultipleAlignment"),
          function(object) {
            return(PairsMAlign$new(object))
          }
)
            
setGeneric("EachWRef", function(object, reference){
  standardGenetic("EachWRef")
})                 
                             
setMethod("EachWRef",
          representation("MultipleAlignment", "character"),
          function(object, reference) {
            return(OneVAll$new(limit = object, one = reference))
          }
)                           
                            




triplets <- function(obj, ...){
  UseMethod('triplets')
}

triplets.DNAMultipleAlignment <- function(obj, checkFunc = function(...) TRUE){
  state <- new.env()
  state$i <- 0L
  state$j <- 1L
  state$k <- 1L
  state$obj <- obj
  it <- list(state=state, checkFunc=checkFunc)
  class(it) <- c("triplets", "containeriter", "iter")
  return(it)
}

getIterVal.triplets <- function(obj, plus = 0L, ...){
  obj$state$i <- obj$state$i + plus
  if(obj$state$i > nrow(obj$state$obj)){
    obj$state$j <- obj$state$j + plus
    obj$state$i <- 0
  }
  if(obj$state$j > nrow(obj$state$obj)){
    obj$state$k <- obj$state$k + plus
    obj$state$j <- 0
  }
  if(obj$state$k > nrow(obj$state$obj) || obj$state$j > nrow(obj$state$obj) || obj$state$i > nrow(obj$state$obj)){
    stop('SubscriptOutOfBounds', call.=FALSE)
  }
  return(obj$state$obj[c(obj$state$i,obj$state$j,obj$state$k)])
}

rangedDataSpaces <- function(obj, split, checkFunc = function(...) TRUE){
  state <- new.env()
  state$i <- 0L
  state$obj <- obj
  it <- list(state = state, checkFunc = checkFunc, recycle = FALSE)
  class(it) <- c("rDataSpaces", "containeriter", "iter", "abstractiter")
  return(it)
}

getIterVal.rDataSpaces <- function(obj, plus = 0L, ...){
  i <- obj$state$i + plus
  if(i > length(obj$state$obj)){stop('SubscriptOutOfBounds', call.=FALSE)}
  return(obj$state$obj[i])
}


# pairsRef <- function(obj, ...){
#   UseMethod('pairsRef')
# }
# 
# pairsRef.DNAMultipleAlignment <- function(obj, ref = NULL, checkFunc = function(...) TRUE){
#   state <- new.env()
#   state$i <- 0L
#   state$obj <- obj
#   if(is.null(ref)){
#     state$ref <- rownames(obj)[1]
#   } else {
#     state$ref <- ref
#   }
#   state$nonRefs <- rownames(obj)
#   state$nonRefs <- state$nonRefs[state$nonRefs != state$ref]
#   it <- list(state=state, checkFunc=checkFunc)
#   class(it) <- c("pairsRef", "containeriter", "iter")
#   return(it)
# }
# 
# getIterVal.pairsRef <- function(obj, plus = 0L, ...){
#   i <- obj$state$i + plus
#   if(i > length(obj$state$nonRefs))
#     stop('SubscriptOutOfBounds', call.=FALSE)
#   pair <- maskSequences(obj$state$obj,
#                         c(obj$state$ref, obj$state$nonRefs[i]),
#                         invert = TRUE,
#                         append = "replace"
#   )
#   return(pair)
# }
# 
# nextElem.pairsRef <- function(obj, ...){
#   repeat {
#     tryCatch({
#       val <- getIterVal(obj, 1L)
#       if(obj$checkFunc(val)){
#         obj$state$i <- obj$state$i + 1L
#         return(val)
#       }
#       obj$state$i <- obj$state$i + 1L
#     }, error = function(e){
#       if(any(nzchar(e$message))){
#         if(identical(e$message, "SubscriptOutOfBounds")){
#           stop("StopIteration", call.=FALSE)
#         } else {
#           stop(e$message, call.=FALSE)
#         }
#       } else {
#         stop("Abort", call.=e)
#       }
#     })
#   }
# }

windows <- function(obj, ...){
  UseMethod('windows')
}

getIterVal <- function(obj, ...){
  UseMethod('getIterVal')
}

windows.default <- function(obj, width, step, inds = FALSE, checkFunc = function(...) TRUE){
  n <- length(obj)
  if(width < 1){stop("Window width must be ≥ 1.")}
  if(step < 1){stop("step must be ≥ 1.")}
  if(any(width > n)){stop("The window size cannot be greater than number of data elements.")}
  
  state <- new.env()
  state$i <- 1L
  state$obj <- obj
  state$width <- width
  state$step <- step
  state$inds <- inds
  
  it <- list(state=state, length=n, checkFunc=checkFunc)
  class(it) <- c("containerwindow", "window", "iter")
  return(it)
}

getIterVal.containerwindow <- function(obj, ...){
  i <- obj$state$i
  if (i > obj$length)
    stop('SubscriptOutOfBounds', call.=FALSE)
  start <- i
  end <- start + obj$state$width - 1
  if(obj$state$inds){
    return(c(start, end))
  } else {
    return(obj$state$obj[start:end])
  }
}

nextElem.containerwindow <- function(obj, ...){
  repeat {
    tryCatch({
      val <- getIterVal(obj)
      if(obj$checkFunc(val)){
        obj$state$i <- obj$state$i + obj$state$step
        return(val)
      }
      obj$state$i <- obj$state$i + obj$state$step
    }, error = function(e){
      if(any(nzchar(e$message))){
        if(identical(e$message, "SubscriptOutOfBounds")){
          stop("StopIteration", call.=FALSE)
        } else {
          stop(e$message, call.=FALSE)
        }
      } else {
        stop("Abort", call.=e)
      }
    })
  }
}