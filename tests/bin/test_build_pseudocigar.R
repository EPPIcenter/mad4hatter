library(testthat)

############################
## Test functions
############################

test_build_pseudoCIGAR_string <- function(build_pseudoCIGAR_string_func) {

  test_that("PseudoCIGAR correctly represents the differences between reference and query", {
    # Example setup
    reference <- "ACT---GT"
    query <- "ACTACTGT"
    expected_output <- "4I=ACT"  # This is just a hypothetical expected output based on your given example

    # Execution
    result <- build_pseudoCIGAR_string_func(reference, query)

    # Assertion
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR outputs a dot when reference and query are the same", {
    # Example setup
    reference <- "ACTGT"
    query <- "ACTGT"
    expected_output <- "."

    # Execution
    result <- build_pseudoCIGAR_string_func(reference, query)

    # Assertion
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR errors out when reference and query are of different lengths", {
    reference <- "ACTG"
    query <- "ACTGT"

    expect_error(build_pseudoCIGAR_string_func(reference, query),
                 "The lengths of reference and query sequences must be the same.")
  })

  test_that("PseudoCIGAR handles single position deletion", {
    reference <- "ACATGT"
    query <- "ACTGT"
    expected_output <- "3D=A"
    expect_error(build_pseudoCIGAR_string_func(reference, query))
  })

  test_that("PseudoCIGAR handles multiple separate insertions and deletions", {
    reference <- "AGCT--GT-CA--T"
    query <- "A-CTAGGTGCAACT"
    expected_output <- "2D=G5I=AG7I=G9I=AC"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR handles multiple separate insertions, masking and deletions", {
    reference <- "AGCT--GT-CA--TGGGTTNNNNNGGGG"
    query <- "A-CTAGGTGCAACT---TTAAAAAG--G"
    expected_output <- "2D=G5I=AG7I=G9I=AC10D=GGG15+5N21D=GG"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR handles masked sequences without changes", {
    reference <- "ACTNNT"
    query <- "ACTGGT"
    expected_output <- "4+2N"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR handles masked sequences with changes", {
    reference <- "ACTNNNNT"
    query <- "ACTG---T"
    expected_output <- "4+4N"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR handles all insertions", {
    reference <- "------"
    query <- "ACTGNT"
    expected_output <- "1I=ACTGNT"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR handles all deletions", {
    reference <- "ACTGNNNNNTA"
    query <- "----------A"
    expected_output <- "1D=ACTG5+5N10D=T"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR handles this monster", {
    reference <- "GGGAC--TGNNNNNTA"
    query <- "GGG---------AA-A"
    expect_warning(build_pseudoCIGAR_string_func(reference, query))
  })

  test_that("PseudoCIGAR handles this monster", {
    reference <- "GGGACTGNNNNNTA"
    query <- "GGG-------AA-A"
    expected_output <- "4D=ACTG8+5N13D=T"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR handles deletions", {
    reference <- "GGGCCCAAA"
    query <- "GGG--C--A"
    expected_output <- "4D=CC7D=AA"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR handles deletions around Ns", {
    reference <- "GGNNNNNGG"
    query <- "G-------G"
    expected_output <- "2D=G3+5N8D=G"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR handles entirely masked reference with indels in middle", {
    reference <- "NN--NN"
    query <- "GGGGGG"
    expected_output <- "1+4N"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })


  test_that("PseudoCIGAR handles entirely masked reference with indels in middle and at end", {
    reference <- "NN--NN--"
    query <- "GGGGGGGG"
    expected_output <- "1+4N"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR handles entirely masked reference with indels in middle and at end with base at end", {
    reference <- "NN--NN--G"
    query <- "GGGGGGGGG"
    expected_output <- "1+4N"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR handles entirely masked reference with indels in middle and at end with insertion at end", {
    reference <- "NN--NN--G-"
    query <- "GGGGGGGGGT"
    expected_output <- "1+4N6I=T"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })


  test_that("PseudoCIGAR handles deletion at first position", {
    reference <- "GGGGGG"
    query <- "-GGGGG"
    expected_output <- "1D=G"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR handles deletion at last position", {
    reference <- "GGGGGG"
    query <- "GGGGG-"
    expected_output <- "6D=G"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR all mask", {
    reference <- "NNNNNN"
    query <- "GGGGG-"
    expected_output <- "1+6N"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })

  test_that("PseudoCIGAR mask on both ends", {
    reference <- "NGGGGGN"
    query <- "GGGGGGG"
    expected_output <- "1+1N7+1N"
    result <- build_pseudoCIGAR_string_func(reference, query)
    expect_equal(result, expected_output)
  })
}