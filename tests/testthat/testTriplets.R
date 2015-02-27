context("tesing functionality of triplets.")

test_that("reference objects storing sequence similarity scan tables work correctly", {
  data(MySequences)
  expect_warning(hyb <- HybRIDS$new(MySequences))
  hyb$analyzeSS(c("Seq1", "Seq2", "Seq3"))
  ssObj <- hyb$triplets$triplets[[1]]$ScanData
  expect_is(ssObj$TableFile, "character")
  expect_is(ssObj$Table, "data.frame")
  expect_true(ssObj$tableIsBlank())
  hyb$analyzeSS(c("Seq1", "Seq2", "Seq3"))
  expect_false(ssObj$tableIsBlank())
  expect_true(!any(apply(ssObj$Table, 2, function(x) any(is.na(x)))))
  ssObj$blankTable()
  expect_true(ssObj$tableIsBlank())
})








