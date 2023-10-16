test_that(".CreateNestedExchangeableCorMat function works",{
  expected_matrix <- matrix(c(1.0, 0.1, 0.1, 0.2, 0.2, 0.2,
                              0.1, 1.0, 0.1, 0.2, 0.2, 0.2,
                              0.1, 0.1, 1.0, 0.2, 0.2, 0.2,
                              0.2, 0.2, 0.2, 1.0, 0.1, 0.1,
                              0.2, 0.2, 0.2, 0.1, 1.0, 0.1,
                              0.2, 0.2, 0.2, 0.1, 0.1, 1.0),
                            nrow = 6, ncol = 6, byrow = TRUE)

  output_matrix <- .CreateNestedExchangeableCorMat(n_obs_periods = 2,
                                                   n_subj_per_period = 3,
                                                   within_period_cor = 0.1,
                                                   between_period_cor = 0.2)

  expect_equal(expected_matrix, output_matrix)
})


# Test case for .CreateExponentialDecayCorMat function
test_that(".CreateExponentialDecayCorMat function works", {
  expected_matrix <- matrix(c(1.000, 0.200, 0.180, 0.180, 0.162, 0.162,
                              0.200, 1.000, 0.180, 0.180, 0.162, 0.162,
                              0.180, 0.180, 1.000, 0.200, 0.180, 0.180,
                              0.180, 0.180, 0.200, 1.000, 0.180, 0.180,
                              0.162, 0.162, 0.180, 0.180, 1.000, 0.200,
                              0.162, 0.162, 0.180, 0.180, 0.200, 1.000),
                            nrow = 6, ncol = 6, byrow = TRUE)

  output_matrix <- round(.CreateExponentialDecayCorMat(n_obs_periods = 3,
                                                       n_subj_per_period = 2,
                                                       within_period_cor = 0.2,
                                                       cor_decay_rate = 0.9), 3)

  expect_equal(expected_matrix, output_matrix)
})

# Test case for .CreateBlockExchangeableCorMat function
test_that(".CreateBlockExchangeableCorMat function works", {
  expected_matrix <- matrix(c(1.0, 0.1, 0.1, 0.3, 0.2, 0.2,
                              0.1, 1.0, 0.1, 0.2, 0.3, 0.2,
                              0.1, 0.1, 1.0, 0.2, 0.2, 0.3,
                              0.3, 0.2, 0.2, 1.0, 0.1, 0.1,
                              0.2, 0.3, 0.2, 0.1, 1.0, 0.1,
                              0.2, 0.2, 0.3, 0.1, 0.1, 1.0),
                            nrow = 6, ncol = 6, byrow = TRUE)

  output_matrix <- .CreateBlockExchangeableCorMat(n_obs_periods = 2,
                                                  n_subj_per_period = 3,
                                                  within_period_cor = 0.1,
                                                  between_period_cor = 0.2,
                                                  within_subject_cor = 0.3)

  expect_equal(expected_matrix, output_matrix)
})

# Test case for .CreateProportionalDecayCorMat function
test_that(".CreateProportionalDecayCorMat function works", {
  expected_matrix <- matrix(c(1.000, 0.200, 0.900, 0.180, 0.810, 0.162,
                              0.200, 1.000, 0.180, 0.900, 0.162, 0.810,
                              0.900, 0.180, 1.000, 0.200, 0.900, 0.180,
                              0.180, 0.900, 0.200, 1.000, 0.180, 0.900,
                              0.810, 0.162, 0.900, 0.180, 1.000, 0.200,
                              0.162, 0.810, 0.180, 0.900, 0.200, 1.000),
                            nrow = 6, ncol = 6, byrow = TRUE)

  output_matrix <- round(.CreateProportionalDecayCorMat(n_obs_periods = 3,
                                                        n_subj_per_period = 2,
                                                        within_period_cor = 0.2,
                                                        cor_decay_rate = 0.9), 3)

  expect_equal(expected_matrix, output_matrix)
})
