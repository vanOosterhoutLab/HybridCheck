# Classes 

# This is a programming construct called an iterator.
# An iterator is some object, typically a small collection of data and some functions, which make it efficient
# to iterate over some container. This makes R's loops more efficient, and better usable in parallel.

eachPair <- function(obj, ...){
  UseMethod('eachPair')
}

eachPair.default <-
  function(obj, checkFunc = function(x, ...) TRUE){
    state <- new.env()
    state$i <- 0L
    state$j <- 0L
    it <- list(state = state, checkFunc = checkFunc)
    class(it) <- c("eachPair", "abstractiter", "iter")
    return(it)
  }

pairsRef <- function(obj, ...){
  UseMethod('pairsRef')
}

pairsRef.default <- 
  function(obj, ref = NULL){
    i <- 0L
    message(sequenceNames(obj))
    if(is.null(ref)){
      ref <- sequenceNames(obj)[1]
    } else {
      ref <- ref
    }
    message(ref, "\n")
    nonRefs <- sequenceNames(obj)
    message(nonRefs)
    nonRefs <- nonRefs[nonRefs != ref]
    message(nonRefs)
    nextEl <- function(){
      i <<- i + 1L
      if(i > length(nonRefs))
        stop('StopIteration', call. = FALSE)
      pair <- subsetSequences(obj,
                              c(ref,
                                nonRefs[i]))
      return(pair)
    }
    it <- list(nextElem = nextEl)
    class(it) <- c("pairsRef", "abstractiter", "iter")
    message("In iterator:")
    message(ref, "\n")
    message(nonRefs, "\n")
    return(it)
  }

#' windows
#' 
#' @export
windows <- function(obj, ...){
  UseMethod('windows')
}

#' windows.default
#' 
#' @export
windows.default <- function(obj, width, step, inds = FALSE, checkFunc = function(x, ...) TRUE){
  n <- length(obj)
  if(width < 1){stop("Window width must be ≥ 1.")}
  if(step < 1){stop("step must be ≥ 1.")}
  if(any(width > n)){stop("The window size cannot be greater than number of data elements.")}
  i <- 1L
  start <- NULL
  end <- NULL
  nextEl <- function(){
    start <<- i
    end <<- start + width - 1
    i <<- i + step
    if (end > length(obj))
      stop('StopIteration', call.=FALSE)
    if(inds){
      return(c(start, end))
    } else {
      return(obj[start:end])
    }
  }
  it <- list(nextElement = nextEl)
  class(it) <- c("containerwindow", "window", "abstractiter", "iter")
  return(it)
}


rangedDataSpaces <- function(obj, split){
  state <- new.env()
  state$i <- 0L
  state$obj <- obj
  it <- list(state = state)
  class(it) <- c("rangedDataSpaces", "abstractiter", "iter")
  return(it)
}

nextElem.rDataSpaces <- function(obj, ...){
  obj$state$i <- obj$state$i + 1L
  if(obj$state$i > length(obj$state$obj)){
    stop('StopIteration', call.=FALSE)
  }
  return(obj$state$obj[i])
}