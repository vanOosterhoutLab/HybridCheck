# Sliding window internal function to zip along sequence triplets and record their sequence similarity as a percentage.
# Created by Ben J. Ward on 18/04/2013.
# Last edited by Ben J. Ward on 02/05/2013.

#' @title Calc.Similarity
#' @description Function performing the sequence similarity analysis with sliding window.
#' @param input An object of type HybRIDSdna.
#' @param window.size A number specifiying the size of the sliding window in base pairs. Must be integer.
#' @param step.size A number specifying the size of the steps the sliding window takes, in base pairs. Must be integer.
#' @export
Calc.Similarity <- function(input, window.size=100, step.size=1){
  if(length(input$Combinations)<2){
    TripletSubset <- list(CroppedSequence=input$CroppedSequence, ContigNames=input$ContigNames)
    analysis <- seq.similarity(TripletSubset, window.size, step.size, fulllength = ncol(input$Sequence), verbose=T)
    return(analysis)
  } else {
    # prepare the outputted list...
    analysis <- list()
    length(analysis) <- length(input$Combinations)
    progress <- txtProgressBar(min=0, max=length(analysis), style=3)
    for(i in 1:length(analysis)){
      setTxtProgressBar(progress, i)
      Triplet <- input$Combinations[[i]]
      TripletSubset <- list(CroppedSequence=input$CroppedSequence[input$Combinations[[i]],], ContigNames=input$ContigNames[Triplet])
      analysis[[i]] <- seq.similarity(TripletSubset, window.size, step.size, fulllength = ncol(input$Sequence), verbose=F)
    }
    return(analysis)
  }
}



seq.similarity <- function(dnain, win.size, s.size, fulllength, verbose) {
  if(verbose==T) cat("\nPreparing input DNA sequences...\n")
  cutDNA <- dnain$CroppedSequence[, colSums(dnain$CroppedSequence[-1,] != dnain$CroppedSequence[-nrow(dnain$CroppedSequence), ]) > 0]
  dnaCutLength <- ncol(cutDNA)     # Measure number of bases in cut data.
  if(dnaCutLength >= 1){   # Only continue with analysis if there is actually some informative sites.
    # First make sure a safe set of sliding windows is calculated for the analysis.
    if(verbose==T) cat("Checking the sliding window parameters\n")
    # 1). Make sure the sliding window size is not some stupid value that will mess up the program loops...
    # Make sure that the window.size value entered by the user is a whole number integer.
    # If not, make it an integer.
    if(is.integer(win.size) == F){          
      win.size <- as.integer(win.size)
    }
    # Make sure the sliding window size is not bigger than the length of the dna sequence.
    # If it is, then automatically set it to 10% of the sequence length.
    # Then if setting the window size to 10% of the sequence size results in an odd window,
    # for example of size 0 or negative correct this by setting the window size to 1.
    # Note: with real data & settings such a arbitrary descision and off setting has never had to be made,
    # by HybRIDS...
    if(win.size > dnaCutLength){
      win.size <- as.integer((dnaCutLength/100)*10)
      if (win.size==0){
        win.size <- 1L
        cat("The set sliding window size is bigger than the length of the actual informative sites of the contig!\n")
        cat("Default behaviour in this case is to set the sliding window to 10%
            of the sequence length, but since this value is below 1, instead HybRIDS is
            setting the sliding window length to 1...\n")
      } else {
        cat("The set sliding window size is bigger than the length of the actual informative sites of the contig!\n")
        cat("Default behaviour in this case is to set the sliding window to 10%
            of the sequence length... \n")
        cat("This is equal to", as.integer((dnaCutLength/100)*10), sep=" ")
      }
    }
    if(verbose==T) cat("Making all the window frames...\n")
    if(win.size >= 1){
      allstepsfrom <- 1 + as.integer(win.size/2)
      allstepsto <- (ncol(cutDNA) - as.integer(win.size/2))+1
      allsteps <- seq(from = allstepsfrom, to = allstepsto, by = s.size)
      windowp1 <- allsteps - as.integer(win.size/2) # All the window start points.
      windowp2 <- allsteps + as.integer(win.size/2) # All the window end points.
      removals <- which(windowp2 > dnaCutLength) # Remove the last window and any accidentally beyond the sequence end point.
      if(length(removals) > 0){
        allsteps <- allsteps[-removals]
        windowp1 <- windowp1[-removals]
        windowp2 <- windowp2[-removals]
      }
      pairs <- combn(1:3,2, simplify=F) # Generate all triplets pairs.
      Distances <- matrix(ncol=9, nrow=length(windowp1))
      Distances[,1] <- allsteps
      Distances[,2] <- windowp1
      Distances[,3] <- windowp2
      Distances[,4] <- as.numeric(unlist(lapply(1:length(allsteps), function(i) colnames(cutDNA)[allsteps[i]])))
      Distances[,5] <- as.numeric(colnames(cutDNA)[windowp1])
      Distances[,6] <- as.numeric(colnames(cutDNA)[windowp2])
      rm(windowp1,windowp2,allsteps,allstepsto,allstepsfrom)
      # Set up the loop for calculation.
      if(verbose==T) cat("Analysing Now!\n")
      if(verbose==T) prog <- txtProgressBar(min=0, max=nrow(Distances), style=3)
      n1 <- 0
      #Do the loop - Calculates all the hammind distances for all contig pairs, in all window frames. 
      for(i in seq(nrow(Distances))) {
        n1 <- n1+1
        dnaStretch <- cutDNA[,Distances[n1,2]:Distances[n1,3]]
        if(verbose==T) setTxtProgressBar(prog, n1)
        Distances[n1,7] <- sum(dnaStretch[pairs[[1]][1],] != dnaStretch[pairs[[1]][2],])
        Distances[n1,8] <- sum(dnaStretch[pairs[[2]][1],] != dnaStretch[pairs[[2]][2],])
        Distances[n1,9] <- sum(dnaStretch[pairs[[3]][1],] != dnaStretch[pairs[[3]][2],])
      }
      colnames( Distances ) <- c( "WindowCenter", "WindowStart", "WindowEnd", "ActualCenter", "ActualStart", "ActualEnd", unlist( lapply( pairs, function( x ) paste( LETTERS[ x ], collapse="" ) ) ) )
      Distances[ , c( 7, 8, 9 ) ] <- 100 - round ( ( as.numeric( Distances[ , c( 7, 8, 9 ) ] ) / ( win.size + 1 ) ) * 100 )
      Distances <- data.frame( Distances )
      } else {
        cat( "The sliding window size is less than 1, this is not supposed to happen...\n" )
        Distances <- "The sliding window size is less than 1, this is not supposed to happen...\n"
        return( Distances )
      }
    } else {
      Distances <- "There are no informative sites to work on - skipping analysis of this triplet...\n"
    }
  return( list(Distances = Distances, InformativeDNAlength = , Contigs = ) )
}


# Internal function For reading in FASTA formatted files.
InputFasta <- function( infile ) {
  dnaSeqLst <- read.fasta( file = infile, 
                           seqtype = "DNA", as.string = FALSE, forceDNAtolower = TRUE,
                           set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = FALSE )
  dnaSeqLengths <- sapply( dnaSeqLst, length )
  if( all( dnaSeqLengths == dnaSeqLengths[[1]] ) == FALSE ) stop( "\nAll sequences are not of the same length. Please check your input files...\nInput files must be fasta files of an alignment, and 'not' the raw sequence files.\nAborting..." )
  cat( "\nAll sequences are of the same length - good - continuing process..." )
  cat( "\nFormatting data..." )
  CompleteDNA <<- do.call( rbind, dnaSeqLst )
  colnames( CompleteDNA ) <<- 1:dnaSeqLengths[[1]]
  rownames( CompleteDNA ) <<- names( dnaSeqLst )
  cat( "\nDNA data formatted successfully!" )
  cat( "\nLooking for duplicates..." )
  dups <- duplicated( CompleteDNA )
  if( any( dups ) ){
    cat("\nSome duplicated sequences were found! - We will get rid of these...")
    CompleteDNA <- CompleteDNA[!dups,]
  }
  return( CompleteDNA )
}