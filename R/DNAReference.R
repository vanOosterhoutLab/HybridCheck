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
                                  FullBp <<- as.numeric( colnames( FullSeq ) )
                                  cat("\nSubsetting the informative segregating sites...")
                                  InformativeSeq <- FullSeq[, colSums( FullSeq[-1,] != FullSeq[-nrow( FullSeq ), ] ) > 0]
                                  InformativeBp <<- as.numeric( colnames( InformativeSeq ) )
                                  SequenceLength <<- ncol( FullSeq )
                                  InformativeLength <<- ncol( InformativeSeq )
                                  SequenceNames <<- rownames( FullSeq )
                                  cat("\nDone, now saving data internally")
                                  cat(" :Full Sequence")
                                  FullSequence <<- FullSeq
                                  cat(" :Informative bases only")
                                  InformativeSequence <<- InformativeSeq
                                  NoDNA <<- FALSE
                                  cat("\nFinished DNA input.")
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
InputSequences <- function( infile, Format ) {
  if(exists(infile)){
    cat("\nFound a variable of that name in the R workspace, now checking it's type...\n")
    classofobj <- class(get(infile))
    if(classofobj == "DNAbin"){
      cat("Detected type of object is DNAbin from the ape package, converting to HybRIDSdna object...\n")
      dna <- as.character(get(infile))
    }
  } else {
    if(is.null(Format)){
      if(!grepl(".", infile)) stop("The provided filename has no extention, you need to specify a format with the forceFormat option.")
      if(grepl(".fas", infile)){
        cat("Detected file is supposed to be a FASTA format file...\n")
        Format <- "fasta"
      }
    }
    cat("Reading in DNA sequence file...")
    dna <- as.character( read.dna( file = infile, format = Format, as.matrix = TRUE ) )
  }
  colnames( dna ) <- 1:ncol( dna )
  cat( "\nLooking for duplicates..." )
  dups <- duplicated( dna )
  if( any( dups ) ){
    cat( "\nSome duplicated sequences were found! - We will get rid of these..." )
    dna <- dna[!dups,]
  }
  cat( "\nDone...")
  return( dna )
}