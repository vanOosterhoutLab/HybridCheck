#' Reference Class for DNA Sequence Information in HybRIDS.
#' @export 
HybRIDSseq <- setRefClass( "HybRIDSseq",
                            fields = list( 
                              FullSequenceFile = "character",
                              FullSequence = function( value ) {
                                if( missing( value ) ){
                                  as.character( read.dna( file = FullSequenceFile, format = "fasta", as.matrix = TRUE ) )
                                } else {
                                  write.dna( value, file = FullSequenceFile, format = "fasta" )
                                }
                              },
                              InformativeSequenceFile = "character",
                              InformativeSequence = function( value ) {
                                if( missing( value ) ){
                                  as.character( read.dna( file = InformativeSequenceFile, format = "fasta", as.matrix = TRUE ) )
                                } else {
                                  write.dna( value, file = InformativeSequenceFile, format = "fasta" )
                                }
                              },
                              SequenceNames = "character",
                              SequenceLength = "numeric",
                              InformativeLength = "numeric",
                              FullBp = "numeric",
                              InformativeBp = "numeric",
                              NoDNA = "logical"),
                              
                            methods = list( 
                              
                              initialize =
                                function( sequenceFile = NULL ) {
                                  FullSequenceFile <<- tempfile( pattern = "FullSequence" )
                                  InformativeSequenceFile <<- tempfile( pattern = "InformativeSequence" )
                                  NoDNA <<- TRUE
                                },
                              
                              InputDNA =
                                function( intarget, format = "FASTA" ) {
                                  if( format == "FASTA" ){
                                    cat("\nLoading in DNA data from FASTA format...")
                                    FullSeq <- InputFasta( intarget )
                                    FullBp <<- as.numeric( colnames( FullSeq ) )
                                    cat("\nMaking finding the informative segregating sites...")
                                    InformativeSeq <- FullSeq[, colSums( FullSeq[-1,] != FullSeq[-nrow( FullSeq ), ] ) > 0]
                                    InformativeBp <<- as.numeric( colnames( InformativeSeq ) )
                                  }
                                  SequenceLength <<- ncol( FullSeq )
                                  InformativeLength <<- ncol( InformativeSeq )
                                  SequenceNames <<- rownames( FullSeq )
                                  cat("\nDone, now saving data internally")
                                  cat(" :Full Sequence")
                                  FullSequence <<- FullSeq
                                  cat(" :Informative bases only")
                                  InformativeSequence <<- InformativeSeq
                                  cat("\nFinished DNA input.")
                                }
                            ))


# Internal function For reading in FASTA formatted files.
InputFasta <- function( infile ) {
  dna <- as.character( read.dna( file = infile, format = "fasta", as.matrix = TRUE ) )
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