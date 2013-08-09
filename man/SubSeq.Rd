\name{SubSeq}
\alias{SubSeq}
\title{Subset a HybRIDSdna object and create a new HybRIDSdna object.}
\usage{
  SubSeq(x, sequences)
}
\arguments{
  \item{x}{A HybRIDSdna object}

  \item{sequences}{Either a vector of sequence numbers, or
  sequence names, the function automaticaly detects which.}
}
\description{
  Create a new HybRIDSdna object from a selection of
  sequences present in another HybRIDSdna object.
}
\details{
  This is a convienience utility function to save the user
  from writing a chunk of code playing the names and
  indexing the dataset. It is also an important function
  for the sequence selector in the HybRIDS GUI.
}

