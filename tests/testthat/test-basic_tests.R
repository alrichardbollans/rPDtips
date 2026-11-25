library(testthat)


test_that("max_pd_tips works correctly", {
  for (i in 1:1000) {
    tr <- rcoal(sample(5:40, 1))
    S  <- sample(tr$tip.label, sample(1:(Ntip(tr) - 1), 1))
    res = max_pd_tips(tr, S)
    expect_true(setequal(res$tip, max_pd_tips_brute(tr, S)$tip))

    phy_diversity <- sum(ape::keep.tip(tr, S)$edge.length)
    gain = sum(ape::keep.tip(tr, c(S, res$tip[1]))$edge.length) - phy_diversity
    expect_equal(round(gain,2), round(res$gain[1],2))

  }

  wrong_input <- "a"
  expect_error(max_pd_tips(wrong_input, NULL))
})







## To test

