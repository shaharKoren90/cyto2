testthat::test_that("Null", {
  testthat::expect_error(
    read_expression_data(),
    "Please provide folder_path or all three seperathe file pathes"
  )

  testthat::expect_error(
    read_expression_data(count_matrix_path='some_path'),
    "Please provide folder_path or all three seperathe file pathes"
  )
})
