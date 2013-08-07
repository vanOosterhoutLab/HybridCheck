\name{dna.data.prepare}
\alias{dna.data.prepare}
\title{A function for loading in DNA sequences from FASTA format and convert to HybRIDS format.}
\usage{
  dna.data.prepare(x = NULL, method = 1, RawThresh = 0.01)
}
\description{
  The function for loading in and converting DNA sequence
  data in FASTA format. It can also perform a pre-analysis
  of whole sequence similarity to try and remove sequences
  unlikely to contain recombinant regions e.g. those
  sequences that may be too similar. Note that analysis
  functions downstream of this function are somewhat
  intelligent in that they know if there are no informative
  sites.
}

