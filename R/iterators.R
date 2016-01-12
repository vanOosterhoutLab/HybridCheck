# Classes 

# This is a programming construct called an iterator.
# An iterator is some object, typically a small collection of data and some functions, which make it efficient
# to iterate over some container. This makes R's loops more efficient, and better usable in parallel.



pairsRef <- function(obj, ...){
  UseMethod('pairsRef')
}

pairsRef.DNAMultipleAlignment <- function(obj, ref = NULL){
  state <- new.env()
  state$i <- 0L
  state$obj <- obj
  if(is.null(ref)){
    state$ref <- rownames(obj)[1]
  } else {
    state$ref <- ref
  }
  state$nonRefs <- rownames(obj)
  state$nonRefs <- state$nonRefs[state$nonRefs != state$ref]
  it <- list(state=state)
  class(it) <- c("pairsRef", "abstractiter", "iter")
  return(it)
}

nextElem.pairsRef <- function(obj, ...){
  obj$state$i <- obj$state$i + 1L
  if(obj$state$i > length(obj$state$nonRefs))
    stop('StopIteration', call.=FALSE)
  pair <- maskSequences(obj$state$obj,
                        c(obj$state$ref, obj$state$nonRefs[obj$state$i]),
                        invert = TRUE,
                        append = "replace"
  )
  return(pair)
}



windows <- function(obj, ...){
  UseMethod('windows')
}

windows.default <- function(obj, width, step, inds = FALSE){
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
  
  it <- list(state=state, length=n)
  class(it) <- c("containerwindow", "window", "abstractiter", "iter")
  return(it)
}

nextElem.containerwindow <- function(obj, ...){
  start <- obj$state$i
  end <- start + obj$state$width - 1
  obj$state$i <- obj$state$i + obj$state$step
  if (end > obj$length)
    stop('StopIteration', call.=FALSE)
  if(obj$state$inds){
    return(c(start, end))
  } else {
    return(obj$state$obj[start:end])
  }
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