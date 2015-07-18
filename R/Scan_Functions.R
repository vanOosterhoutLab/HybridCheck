# Functions for sliding window analyses.

#' @title Check sliding window settings used for a scan.
#' @description Internal function. Used to make sure the user has entered
#' a sensible option for the sliding window size and step increments.
#' @param winSize A single numeric value. The proposed size of the sliding
#' window in the scan.
#' @param trackLen The A single numeric value. length of the track that the
#' window will travel down.
#' @examples 
#' HybridCheck:::windowSizeChecker(10, 100)
#' HybridCheck:::windowSizeChecker(100, 50)
#' HybridCheck:::windowSizeChecker(100, 10)
#' HybridCheck:::windowSizeChecker(100, 3)
#' @return winSize Single numeric value. The size of the sliding window to
#' be used in the analysis after the check. If an accepable size was chosen by
#' the user, this value will be identical to the input value. 
windowSizeChecker <- function(winSize, trackLen){
  message("\t- Checking the sliding window parameters.")
  if(winSize > trackLen){
    winSize <- as.integer((trackLen / 100) * 10)
    warning("\t\t- The set sliding window size is bigger than the length of the actual informative sites of the contig!")
    message("\t\t- Continuing with analysis but set the sliding window to 10%
            of the sequence length... ")
    message("\t\t- This is equal to ", winSize)
    if(winSize < 1L){
      winSize <- 1L
      message("\t\t- Set the sliding window to 10% of the sequence length, 
              but since this value is below 1, instead setting
              the sliding window length to 1...")
      }
    }
  return(winSize)
}

#' @title Generate data-frame with co-ordinates of windows for analysis.
#' @description Internal function. Used to make a data-frame of all the sliding
#' windows: their mid-points, start, ends. This data-frame is used during the 
#' sequence alignment scans and also forms part of the final data table that
#' is produced.
#' @param winSize Single numeric value. The size (in bp) of the sliding window.
#' @param stepSize Single numeric value. The number of base positions the 
#' sliding window jumps along the sequences every time it moves.
#' @param trackLen Single numeric value. Length of the track the window will 
#' travel down.
#' @param bases A numeric vector. The base positions that the sliding window
#' will slide over. Note in many HybridCheck analyses, the sliding window only
#' slides over the informative sites in a sequence alignment i.e. fully
#' conserved, non-informative sites are removed.
#' @return slideFrame A data-frame of six columns. Each column is numeric.
makeWindowFrames <- function(winSize, stepSize, trackLen, bases){
  message("\t- Making all the window frames...")
  if(winSize >= 1L){
    halfWindow <- as.integer(winSize / 2)
    allstepsfrom <- 1 + halfWindow
    allstepsto <- (trackLen - halfWindow) + 1
    allsteps <- seq(from = allstepsfrom, to = allstepsto, by = stepSize)
    windowp1 <- allsteps - halfWindow # All the window start points.
    windowp2 <- allsteps + halfWindow # All the window end points.
    removals <- which(windowp2 > trackLen)
    if(length(removals) > 0) {
      allsteps <- allsteps[-removals]
      windowp1 <- windowp1[-removals]
      windowp2 <- windowp2[-removals]
    }
    slideFrame <- matrix(ncol = 6, nrow = length(windowp1))
    slideFrame[, 1] <- allsteps
    slideFrame[, 2] <- windowp1
    slideFrame[, 3] <- windowp2
    # ActualBP Center
    slideFrame[, 4] <- as.numeric(unlist(lapply(1:length(allsteps), function(i) bases[allsteps[i]])))
    # Actual BP Start
    slideFrame[, 5] <- as.numeric(bases[windowp1])
    # Actual BP End
    slideFrame[, 6] <- as.numeric(bases[windowp2]) 
    # Make placeholders for scan results.
    return(slideFrame)
  } else {
    stop("The sliding window size is less than 1, this is not supposed to be possible.")
  }
}

makeConMats <- function(dnaSequences, pairs){
  conMats <- lapply(pairs, function(i){
    colSums(consensusMatrix(dnaSequences[i]) != 0) > 1
  })
  return(conMats)
}

calculateDistanceTracks <- function(dnaSequences, pairs, resTable){
  consensusMatrices <- makeConMats(dnaSequences, pairs)
  mutationCountTracks <- lapply(consensusMatrices, function(x){
    unlist(lapply(seq(nrow(resTable)), function(i){
      stretch <- resTable[i, 2] : resTable[i, 3]
      return(sum(x[stretch]))
    }))
  })
  return(mutationCountTracks)
}


scan.similarity <- function(dna, triplet, ambiguousAreHet, settings){
  message(paste0(" - Scanning sequence similarity for triplet ",
                 paste0(triplet$SequenceInfo$ContigNames, collapse=", ")))
  cutDNA <- triplet$SequenceInfo$prepareDNAForScan(dna, ambiguousAreHet)
  triplet$readSettings(settings)
  if(triplet$SequenceInfo$InformativeUsedLength >= 1){
    triplet$ScanData$WindowSizeUsed <- 
      windowSizeChecker(triplet$ScanData$WindowSizeUsed,
                        triplet$SequenceInfo$InformativeUsedLength)
    if(triplet$ScanData$WindowSizeUsed >= 1L) {
      sequencePairs <- combn(1:3 , 2, simplify = F)
      Distances <- makeWindowFrames(triplet$ScanData$WindowSizeUsed,
                                    triplet$ScanData$StepSizeUsed,
                                    triplet$SequenceInfo$InformativeUsedLength,
                                    triplet$SequenceInfo$InformativeUsed,
                                    sequencePairs)
      # Set up the loop for calculation.
      message("\t- Scanning Now!")
      tracks <- calculateDistanceTracks(cutDNA, sequencePairs, Distances)
      Distances <- cbind(Distances, do.call(cbind, tracks))
      colnames(Distances) <- c("WindowCenter", "WindowStart", "WindowEnd", 
                               "ActualCenter", "ActualStart", "ActualEnd",
                               unlist(
                                 lapply(sequencePairs, 
                                        function(x) paste(LETTERS[x], 
                                                          collapse=""))))
      Distances[ , c(7, 8, 9)] <- 100 - round((as.numeric(Distances[ , c(7, 8, 9)]) / (triplet$ScanData$WindowSizeUsed + 1)) * 100)
      triplet$ScanData$Table <- as.data.frame(Distances)
    } else {
      stop("The sliding window size is less than 1, this is not supposed to be possible.")
    }
  } else {
    warning(paste0("There are no informative sites to work on - skipping analysis of this triplet: ", triplet$ContigNames, collapse = ", "))
  }
}

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