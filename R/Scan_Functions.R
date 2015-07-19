# Functions for sliding window analyses.






# Get the bases which the heterozygous codes fed in, have in common.
inCommon <- function(code, bases){
  return(Reduce(intersect, code[bases]))
}

# Get the heterozygous sites located at a particular base position.
getAmbig <- function(mat, sel){
  return(names(which(mat[5:10, sel] > 0)))
}

transformSequence <- function(dnaSeq, transTable){
  sitesToBeTrans <- dnaSeq[as.numeric(transTable$Base)]
  matchesToTransTable <- cbind(
    strsplit(as.character(sitesToBeTrans), "")[[1]] == transTable[,2],
    strsplit(as.character(sitesToBeTrans), "")[[1]] == transTable[,3],
    strsplit(as.character(sitesToBeTrans), "")[[1]] == transTable[,4])
  matchesToTransTable[is.na(matchesToTransTable)] <- FALSE
  transTo <- numeric(length=nrow(matchesToTransTable))
  transTo[matchesToTransTable[,1]] <- 1
  transTo[matchesToTransTable[,2]] <- 2
  transTo[matchesToTransTable[,3]] <- 3
  transTo[transTo == 0] <- NA
  sitesToTrans <- transTable$Base[which(!is.na(transTo))]
  transTable <- as.matrix(transTable[which(!is.na(transTo)),5:7])
  midx <- cbind(1:nrow(transTable), transTo[which(!is.na(transTo))])
  transSeq <- paste0(transTable[midx], collapse = "")
  return(replaceLetterAt(dnaSeq, sitesToTrans, transSeq))
}

transformSequenceRelative <- function(dnaSeq, transTable, firstBase, lastBase){
  transInSegment <- 
    transTable[which(
      (transTable$TrueBase >= firstBase) &
        (transTable$TrueBase <= lastBase)),]
  transInSegment$Base <- (transInSegment$TrueBase - firstBase) + 1
  transformSequence(dnaSeq, transInSegment)
}


# Identify and decide on transformations for sites with one type of heterozygous base.
transSingleAmb <- function(code, atypes, statemat){
  oneBase <- which(atypes == 1)
  if(length(oneBase) > 0){
    # - Figure out what the ambiguous state is in each case...
    oneLetters <- unlist(lapply(oneBase, function(x) getAmbig(statemat, x)))
    # - Decide on what bases the ambigious sites should be transformed to...
    oneTrans <- lapply(oneLetters, function(i) sample(code[[i]], 1))
    return(data.frame(Base = unlist(oneBase), AmbigOne = unlist(oneLetters), AmbigTwo = NA,
                      AmbigThree = NA, ResolveOne = unlist(oneTrans), ResolveTwo = NA,
                      ResolveThree = NA))
  } else {
    return(data.frame(Base = NA, AmbigOne = NA, AmbigTwo = NA, AmbigThree = NA,
                      ResolveOne = NA, ResolveTwo = NA, ResolveThree = NA))
  }
}

# Identify and decide on transformations of sites with two kinds of hetrozygous base.
transTwoAmb <- function(code, atypes, statemat){
  # Figure out which bases are bases with two kinds of ambiguity.
  twoBase <- which(atypes == 2)
  if(length(twoBase) > 0){
    # Get the two ambiguous bases at each site identified.
    twoLetters <- lapply(twoBase, function(x){getAmbig(statemat, x)})
    # Decide randomly which of the two ambiguous sites to pick to resolve first.
    whichOfTwo <- sample(c(1, 2), length(twoLetters), replace = T)
    chosenAmbig <- mapply(function(x, y){
      x[y]
    }, x = twoLetters, y = whichOfTwo)
    otherAmbig <- mapply(function(x, y){
      x[which(x != y)]
    }, twoLetters, chosenAmbig)
    # For each chosen ambiguous base, randomly pick one of its bases from the code.
    chosenBases <- unlist(lapply(chosenAmbig, function(x){sample(code[[x]], 1)}))
    # For the second ambiguous base, if the chosen base picked previously is also encoded for
    # by the ambiguous base, then resolve it the same way, otherwise, pick from its code randomly.
    chosenBases2 <- mapply(function(ambigs, pickedbase, altambig){
      if(pickedbase %in% inCommon(code, ambigs)){
        return(pickedbase)
      } else {
        sample(code[[altambig]], 1)
      }
    }, twoLetters, chosenBases, otherAmbig)
    # Collect the results into a datastructure and return it.
    return(data.frame(Base = twoBase, AmbigOne = chosenAmbig, AmbigTwo = otherAmbig,
                      AmbigThree = NA, ResolveOne = chosenBases, 
                      ResolveTwo = chosenBases2, ResolveThree = NA))
  } else {
    return(data.frame(Base = NA, AmbigOne = NA, AmbigTwo = NA, AmbigThree = NA,
                      ResolveOne = NA, ResolveTwo = NA, ResolveThree = NA))
  }
}

# Identify and decide on transformations for sites with three distinct kinds of heteroztgous base.
transThreeAmb <- function(code, atypes, statemat){
  threeBase <- which(atypes == 3)
  if(length(threeBase) > 0){
    # Get the three ambiguous bases at each site identified.
    threeLetters <- lapply(threeBase, function(x){getAmbig(statemat, x)})
    oneToDiscard <- lapply(threeLetters, function(x){
      return(sample(x, 1))
    })
    leftToResolve <- mapply(function(x, y){
      return(x[x != y])
    }, threeLetters, oneToDiscard, SIMPLIFY = FALSE)
    resolveDiscard <- lapply(leftToResolve, function(x){
      return(sample(x, 1))
    })
    twoLetters <- mapply(function(x, y, z){
      x[which(x == y)] <- z
      return(unique(x))
    }, threeLetters,  oneToDiscard, resolveDiscard, SIMPLIFY = FALSE)
    whichOfTwo <- sample(c(1, 2), length(twoLetters), replace = T)
    chosenAmbig <- mapply(function(x, y){
      x[y]
    }, x = twoLetters, y = whichOfTwo)
    otherAmbig <- mapply(function(x, y){
      x[which(x != y)]
    }, twoLetters, chosenAmbig)
    chosenBases <- lapply(chosenAmbig, function(x){sample(code[[x]], 1)})
    chosenBases2 <- mapply(function(ambigs, pickedbase, altambig){
      if(pickedbase %in% inCommon(code, ambigs)){
        return(pickedbase)
      } else {
        sample(code[[altambig]], 1)
      }
    }, twoLetters, chosenBases, otherAmbig)
    return(data.frame(Base = unlist(threeBase), AmbigOne = unlist(chosenAmbig), AmbigTwo = unlist(otherAmbig),
                      AmbigThree = unlist(oneToDiscard), ResolveOne = unlist(chosenBases),
                      ResolveTwo = unlist(chosenBases2), ResolveThree = unlist(resolveDiscard)))
  } else {
    return(data.frame(Base = NA, AmbigOne = NA, AmbigTwo = NA, AmbigThree = NA,
                      ResolveOne = NA, ResolveTwo = NA, ResolveThree = NA))
  }
}