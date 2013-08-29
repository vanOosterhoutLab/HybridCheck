# Reference Class for DNA Sequence Information in HybRIDS 
DNASequence <- setRefClass( "HybRIDSseq",
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
                              InformativeBp = "numeric" ),
                              
                            methods = list( 
                              
                              initialize =
                                function( sequenceFile = NULL, fileFormat = "FASTA") {
                                  if( !is.null( sequenceFile ) ) {
                                    InputDNA( sequenceFile, fileFormat )
                                  }
                                },
                              
                              InputDNA =
                                function( intarget, format = "FASTA" ) {
                                  FullSequenceFile <<- tempfile( pattern = "FullSequence" )
                                  InformativeSequenceFile <<- tempfile( pattern = "InformativeSequence" )
                                  if( format == "FASTA" ){
                                    FullSeq <- InputFasta( intarget )
                                    FullBp <<- as.numeric( colnames( FullSeq ) )
                                    InformativeSeq <- FullSeq[, colSums( FullSeq[-1,] != FullSeq[-nrow( FullSeq ), ] ) > 0]
                                    InformativeBp <<- as.numeric( colnames( InformativeSeq ) )
                                    FullSequence <<- FullSeq
                                    InformativeSequence <<- InformativeSeq
                                  }
                                  SequenceNames <<- rownames( FullSequence )
                                  SequenceLength <<- ncol( FullSequence )
                                },
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
  return( dna )
}