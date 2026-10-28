skip_unless_slow_tests <- function() {
  testthat::skip_if(
    !identical(
      Sys.getenv("BFBIN2ARM_RUN_SLOW_TESTS"),
      "true"
    ),
    paste(
      "Slow numerical validation test.",
      "Set BFBIN2ARM_RUN_SLOW_TESTS=true to run."
    )
  )
}