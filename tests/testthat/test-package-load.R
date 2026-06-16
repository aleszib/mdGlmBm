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
    "summary.dynamic_network",
    "fit_time_glm_blockmodel",
    "fit_time_glm_blockmodels",
    "estimate_markov_transitions",
    "estimate_membership_prior",
    "score_actor_time_candidates",
    "fit_dynamic_glm_blockmodel",
    "print.time_glm_blockmodel",
    "summary.time_glm_blockmodel",
    "print.time_glm_blockmodels",
    "summary.time_glm_blockmodels",
    "print.markov_transitions",
    "summary.markov_transitions",
    "print.membership_prior",
    "summary.membership_prior",
    "print.dynamic_glm_blockmodel",
    "summary.dynamic_glm_blockmodel"
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
