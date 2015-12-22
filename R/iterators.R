# Classes 

# This is a programming construct called an iterator.
# An iterator is some object, typically a small collection of data and some functions, which make it efficient
# to iterate over some container. This makes R's loops more efficient, and better usable in parallel.


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


pairsRef <- function(obj, ...){
  UseMethod('pairsRef')
}

pairsRef.DNAMultipleAlignment <- function(obj, ref = NULL, checkFunc = function(...) TRUE){
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
  it <- list(state=state, checkFunc=checkFunc)
  class(it) <- c("pairsRef", "containeriter", "iter")
  return(it)
}

getIterVal.pairsRef <- function(obj, plus = 0L, ...){
  i <- obj$state$i + plus
  if(i > length(obj$state$nonRefs))
    stop('SubscriptOutOfBounds', call.=FALSE)
  pair <- maskSequences(obj$state$obj,
                        c(obj$state$ref, obj$state$nonRefs[i]),
                        invert = TRUE,
                        append = "replace"
  )
  return(pair)
}

nextElem.pairsRef <- function(obj, ...){
  repeat {
    tryCatch({
      val <- getIterVal(obj, 1L)
      if(obj$checkFunc(val)){
        obj$state$i <- obj$state$i + 1L
        return(val)
      }
      obj$state$i <- obj$state$i + 1L
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