# #' @include Generic_Methods.R
# NULL
# 
# setMethod("sequenceLength",
#           representation("DNAMultipleAlignment"),
#           function(object) {
#             return(ncol(object))
#           })
# 
# setMethod("maskConservedSites",
#           signature("DNAMultipleAlignment", "character"),
#           function(object, append) {
#             conserved <-
#               IRanges(na.omit(colSums(consensusMatrix(object) != 0) == 1))
#             colmask(object, append = append) <- conserved
#             return(object)
#           })
# 
# setMethod("excludeUnknownSites",
#           signature("DNAMultipleAlignment", "character"),
#           function(object, append) {
#             ns <- IRanges(consensusMatrix(object)[15, ] != 0)
#             colmask(object, append = append) <- ns
#             return(object)
#           })
# 
# setMethod("maskedSites",
#           signature("DNAMultipleAlignment"),
#           function(object) {
#             return(which(as.logical(coverage(colmask(
#               object
#             )))))
#           })
# 
# setMethod("unmaskedSites",
#           signature("DNAMultipleAlignment"),
#           function(object) {
#             (1:sequenceLength(object))[-maskedSites(object)]
#           })
# 
# setMethod("nMaskedSites",
#           signature("DNAMultipleAlignment"),
#           function(object) {
#             sum(colmask(object)@width)
#           })
# 
# setMethod("nUnMaskedSites",
#           signature("DNAMultipleAlignment"),
#           function(object) {
#             ncol(object) - nMaskedSites(object)
#           })
# 
# setMethod("subsetSequences",
#           signature("DNAMultipleAlignment", "character", "logical", "character"),
#           function(object, seqnames, invert, append) {
#             sequenceNames <- rownames(object)
#             rowmask(object, append = append, invert = invert) <-
#               which(sequenceNames %in% seqnames)
#             return(object)
#           })
# 
# setMethod("getSeqNames",
#           signature("DNAMultipleAlignment"),
#           function(object) {
#             rownames(object)[-as.integer(rowmask(object))]
#           })
# 
# 
# setMethod("statesPerBase",
#           signature("DNAMultipleAlignment"),
#           function(object) {
#             return(colSums(consensusMatrix(object) != 0))
#           })
# 
# setMethod("polySites",
#           signature("DNAMultipleAlignment"),
#           function(object) {
#             return(statesPerBase(object) > 1)
#           })
# 
# setMethod("slidingPoly",
#           signature("DNAMultipleAlignment", "integer", "integer"),
#           function(object, windowSize, stepSize) {
#             conMat <- polySites(object)
#             itr <- windows(
#               conMat,
#               width = windowSize,
#               step = stepSize,
#               checkFunc = function(i)
#                 ! any(is.na(i))
#             )
#             dists <- foreach(x = itr, .combine = c) %dopar% {
#               sum(x)
#             }
#             dists <- 100 - round((dists / windowSize) * 100)
#             windowStarts <-
#               seq(from = 1, to = nUnMaskedSites(object), by = stepSize)
#             windowEnds <-
#               seq(from = windowSize, to = nUnMaskedSites(object), by = stepSize)
#             data <-
#               RangedData(
#                 space = paste(getSeqNames(object), collapse = ":"),
#                 ranges = IRanges(start = windowStarts[1:length(windowEnds)],
#                                  end = windowEnds),
#                 signal = dists
#               )
#             return(data)
#           })