test_that("truncation occurs coordinatewise before covariance construction", {
  block <- matrix(c(
    1, 2, 3,
    4, 5, 6,
    7, 8, 9
  ), nrow = 3, byrow = TRUE)
  thresholds <- c(5, 6, 8)
  truncated <- hfsis:::.truncate_coordinates(block, thresholds)
  expected <- block
  expected[abs(block) > matrix(thresholds, 3, 3, byrow = TRUE)] <- 0
  expect_equal(truncated, expected)

  covariance <- hfsis:::.local_truncated_covariance(
    block, target_index = 1L, k = 3L, delta = 0.5,
    thresholds = thresholds
  )
  expect_equal(covariance, crossprod(expected) / (3 * 0.5))
  expect_gte(min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values),
    -1e-12)
})
