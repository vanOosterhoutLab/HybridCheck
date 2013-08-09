\name{SimplifyBlocks}
\alias{SimplifyBlocks}
\alias{SimplifyBlocks.HybRIDSblock}
\alias{SimplifyBlocks.HybRIDSblockSET}
\alias{SimplifyBlocks.HybRIDSdatedBlocks}
\alias{SimplifyBlocks.HybRIDSdatedBlocksSET}
\title{SimplifyBlocks}
\usage{
  SimplifyBlocks(x, ...)

  SimplifyBlocks.HybRIDSblock(x)

  SimplifyBlocks.HybRIDSdatedBlocks(x)

  SimplifyBlocks.HybRIDSblockSET(x, seq1,
    seq2 = "Individual", seq3 = "Pair", returnEach = T,
    selectionOnly = T)

  SimplifyBlocks.HybRIDSdatedBlocksSET(x, seq1,
    seq2 = "Individual", seq3 = "Pair", returnEach = T,
    selectionOnly = T)
}
\arguments{
  \item{x}{an object of class HybRIDSblock,
  HybRIDSdatedBlock, HybRIDSblockSET, or
  HybRIDSdatedBlocksSET.}
}
\description{
  A function to simplify HybRIDSblock, HybRIDSdatedBlock,
  HybRIDSblockSET, and HybRIDSdatedBlocksSET objects as
  tables which can be exported.
}

