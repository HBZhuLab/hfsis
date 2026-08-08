test_that("example data use the documented response key", {
  data(hf_sis_example, package = "hfsis")
  expect_identical(hf_sis_example$response, "y")
  expect_true(hf_sis_example$response %in% names(hf_sis_example$data))
})
