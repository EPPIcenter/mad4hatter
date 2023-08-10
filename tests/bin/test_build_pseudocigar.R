library(testthat)

module <- file.path(Sys.getenv("TESTDIR"), "build_pseudocigar.R")
source(module) # source the module to test

test_that("build_pseudoCIGAR_string with UNMASKED data", {
  reference <- "AAGCTTTGCA-"
  query <- "AAGCTTGCA-"
  expect_equal(build_pseudoCIGAR_string(reference, query), "5T6I=T8D=C")
})

test_that("build_pseudoCIGAR_string with MASKED data", {
  reference <- "AAGCTTTGCA-N"
  query <- "AAGCTTGCA-N"
  expect_equal(build_pseudoCIGAR_string(reference, query), "5T6I=T8D=C")
})

