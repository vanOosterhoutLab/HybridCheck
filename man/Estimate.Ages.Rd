\name{Estimate.Ages}
\alias{Estimate.Ages}
\title{Function that estimates the ages of putative recombination blocks}
\usage{
  Estimate.Ages(blocks, sequence, mutation.rate = 1e-07,
    requiredP = 0.005)
}
\arguments{
  \item{blocks}{An object of class HybRIDSblock or
  HybRIDSblockSET}

  \item{sequence}{The object of class HybRIDSdna that was
  used to generate the object given as the 'blocks'
  parameter.}
}
\description{
  This function take a HybRIDSblock or HybRIDSblockSET and
  will estimate the ages of the blocks detected.
}

