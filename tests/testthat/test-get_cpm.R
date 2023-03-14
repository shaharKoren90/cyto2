testthat::test_that("Wrong formay", {
  testthat::expect_error(
    get_cpm('a_string'),
    "Data should be a count matrix or an ExpressionSet"
  )

  testthat::expect_error(
    get_cpm('a_vector'),
    "Data should be a count matrix or an ExpressionSet"
  )
})
