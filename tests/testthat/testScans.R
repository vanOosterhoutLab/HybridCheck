context("Sliding window scans")

test_that("The function windowSizeChecker correctly evaluates and modifies 
          Sliding window size values chosen by the user.", {
            expect_equal(HybridCheck:::windowSizeChecker(10, 100), 10)
            expect_message(HybridCheck:::windowSizeChecker(10, 100))
            expect_equal(HybridCheck:::windowSizeChecker(100, 50), 5)
            expect_warning(HybridCheck:::windowSizeChecker(100, 50))
            expect_equal(HybridCheck:::windowSizeChecker(100, 10), 1)
            expect_warning(HybridCheck:::windowSizeChecker(100, 10))
            expect_equal(HybridCheck:::windowSizeChecker(100, 3), 1)
            expect_warning(HybridCheck:::windowSizeChecker(100, 3))
            })

test_that("The function makeWindowFrames corretly creates a data-frame of
          slidingWindow co-ordinates and make sure contents of the frame
          make sense.", {
            bases <- 1:1000
            windowSize <- stepSize <- 100
            trackLength <- 1000
            testFrame <- HybridCheck:::makeWindowFrames(windowSize, stepSize, trackLength, bases)
            expect_false(nrow(testFrame) > trackLength)
            expect_true(all(
              testFrame[, 1] > testFrame[, 2] &
                testFrame[, 1] < testFrame[, 3] &
                testFrame[, 3] > testFrame[, 1]))
            expect_true(all(
              testFrame[, 4] > testFrame[, 5] &
                testFrame[, 4] < testFrame[, 6] &
                testFrame[, 6] > testFrame[, 4]))
            expect_error(HybridCheck:::makeWindowFrames(-2, stepSize, trackLength, bases))
          })
