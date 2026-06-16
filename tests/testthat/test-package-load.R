test_that("package loads and exported functions exist", {
  expect_silent(library("dynGLMbm", character.only = TRUE))

  exports <- c(
    "mdsbm_opr",
    "mdsbm_one_partition",
    "mdsbm_icl_one_partition",
    "ppml",
    "optParGlm",
    "optRandomParGlm",
    "optRandomParRangeGlm",
    "upAndDownSearch",
    "balassaNorm",
    "glm_blockmodel_family",
    "as_dynamic_network",
    "print.dynamic_network",
    "summary.dynamic_network"
  )

  for (fn in exports) {
    expect_true(exists(fn, envir = asNamespace("dynGLMbm"), inherits = FALSE))
  }
})

test_that("no package R files contain top-level library calls", {
  r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
  for (file in r_files) {
    lines <- readLines(file, warn = FALSE)
    expect_false(any(grepl("^\\s*library\\(", lines)), info = file)
  }
})
