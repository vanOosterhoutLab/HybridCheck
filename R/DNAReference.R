#' Reference Class for DNA Sequence Information in HybRIDS.
#' @export 
HybRIDSseq <- setRefClass( "HybRIDSseq",
                            fields = list( 
                              FullSequence = "ANY",
                              InformativeSequence = "ANY",
                              SequenceNames = "character",
                              SequenceLength = "numeric",
                              InformativeLength = "numeric",
                              FullBp = "numeric",
                              InformativeBp = "numeric",
                              InformativeLoci = "numeric",
                              NoDNA = "logical"),
                              
                            methods = list( 
                              initialize =
                                function(sequenceInput = NULL, force_Format = NULL) {
                                  NoDNA <<- TRUE
                                  if(!is.null(sequenceInput)){
                                    InputDNA(sequenceInput, force_Format)
                                  }
                                },
                              
                              InputDNA =
                                function( intarget, forceFormat = NULL) {
                                  FullSeq <- InputSequences(intarget, forceFormat)
                                  FullBp <<- as.numeric(colnames(FullSeq))
                                  message("Subsetting the informative segregating sites...")
                                  InformativeSeq <- FullSeq[, colSums(FullSeq[-1,] != FullSeq[-nrow(FullSeq), ] ) > 0]
                                  InformativeBp <<- as.numeric(colnames(InformativeSeq))
                                  SequenceLength <<- ncol(FullSeq)
                                  InformativeLength <<- ncol(InformativeSeq)
                                  SequenceNames <<- rownames(FullSeq)
                                  message("Done, now saving data internally")
                                  message(" :Full Sequence")
                                  FullSequence <<- FullSeq
                                  message(" :Informative bases only")
                                  InformativeSequence <<- InformativeSeq
                                  NoDNA <<- FALSE
                                  message("Finished DNA input.")
                                },
                              
                              plotInf =
                                function(parameters, which="hist") {
                                  if(which != "hist" && which != "bars") stop("You need to specify either 'bars' or 'hist' as the which option!")
                                  if (which == "bars"){
                                    # Figure out the scale and data to go into each vertical bar: TODO - Put this in a function.
                                    div <- FullBp[SequenceLength] / parameters$MosaicScale
                                    bpstart <- seq(from = 1, to = FullBp[SequenceLength], by = div)
                                    bpend <- seq(from=div, to = FullBp[SequenceLength], by = div)
                                    bpX <- round(bpstart + (div / 2))
                                    frame <- data.frame(bpstart = bpstart, bpend = bpend, bpX = bpX)
                                    rm(bpstart, bpend)
                                    Bar <- round(apply(frame, 1, function(x) vertbar_create2(InformativeBp, x)))
                                    plotFrame <- data.frame(bpstart = frame$bpstart, bpend = frame$bpend, X = frame$bpX, numberSNPs = Bar)
                                    thePlot <- ggplot(plotFrame, aes(x = X, y = 1)) + geom_raster(aes(fill = numberSNPs))
                                    thePlot <- applyPlottingParams(thePlot, parameters, title="Number of SNP Polymorphisms Between All Sequences")
                                    thePlot <- thePlot + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
                                  } else {
                                    thePlot <- ggplot(data.frame(InformativeBP = InformativeBp), aes(x = InformativeBP)) + geom_histogram(aes(fill=..count..)) +
                                      scale_fill_gradient2("SNP Count", low="white", high="red")
                                    thePlot <- applyPlottingParams(thePlot, parameters, title="Number of SNP Polymorphisms Between All Sequences")
                                    thePlot <- thePlot + ylab("Count")
                                  }
                                  thePlot <- thePlot + xlab("Base Pair Position")
                                  return(thePlot)
                                }
                            ))


HybRIDSseq_fastadisk <- setRefClass( "HybRIDSseq_fastadisk",
                                     contains = "HybRIDSseq",
                           fields = list( 
                             FullSequenceFile = "character",
                             FullSequence = function( value ) {
                               if( missing( value ) ){
                                 as.character( read.dna( file = FullSequenceFile, format = "fasta", as.matrix = TRUE ) )
                               } else {
                                 write.dna( as.DNAbin(value), file = FullSequenceFile, format = "fasta" )
                               }
                             },
                             InformativeSequenceFile = "character",
                             InformativeSequence = function( value ) {
                               if( missing( value ) ){
                                 as.character( read.dna( file = InformativeSequenceFile, format = "fasta", as.matrix = TRUE ) )
                               } else {
                                 write.dna( value, file = InformativeSequenceFile, format = "fasta" )
                               }
                             }),
                           
                           methods = list( 
                             initialize =
                               function(sequenceInput = NULL, force_Format = NULL){
                                 FullSequenceFile <<- tempfile( pattern = "FullSequence" )
                                 InformativeSequenceFile <<- tempfile( pattern = "InformativeSequence" )
                                 if(!is.null(sequenceInput)){
                                   InputDNA(sequenceInput, force_Format)
                                 }
                               }
                           ))


# Internal function For reading in sequence files, based on the format deteted.
InputSequences <- function(input, Format) {
  classOfInput <- class(input)
  if(classOfInput == "character"){
    # Class of input is text, so we assume it is a filepath...
    # We need to check the format of the file to be read in and then indeed read it in.
    if(grepl(".fas", input) || Format == "FASTA" || Format == "fasta"){
      message("File to be read is expected to be FASTA format...")
      Format <- "fasta"
    }
    message("Reading in sequence file...")
    dna <- read.dna( file = input, format = Format, as.matrix = TRUE )
  } else {
    if(classOfInput == "DNAbin"){
      message("Class of input is DNAbin from ape package.")
      dna <- input
    }
  }
  message("Looking for duplicates...")
  distances <- dist.dna(dna, model = "N")
  if(any(distances == 0)){
    message("Some duplicated sequences were found! - We will get rid of these...")
    indicies <- distrowcol(which(distances == 0), attr(distances, "Size"))
    dna <- dna[-indicies[,2],]
    message("Double Checking for duplicated sequences again to be safe...")
    distances <- dist.dna(dna, model="N")
    if(any(distances == 0)) stop("Duplicates were still found in the sequences - this should not happen - aborting. Inform package maintainer.")
  }
  dna <- as.character(dna)
  colnames(dna) <- 1:ncol(dna)
  message("Done...")
  return(dna)
}
  
    
distrowcol <- function(ix,n){
  nr <- ceiling(n-(1+sqrt(1+4*(n^2-n-2*ix)))/2)
  nc <- n-(2*n-nr+1)*nr/2+ix+nr
  return(cbind(nr,nc))
}