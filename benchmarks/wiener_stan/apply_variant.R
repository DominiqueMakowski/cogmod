# Patch a working tree so that the DDM-4 and DDM-5 call site of
# cogmod_ddm_decision_lpdf() uses a Stan-language Wiener density from this
# directory instead of Stan's wiener_lpdf(). The 7-parameter call site, the
# router it uses and everything on the R side are left alone. Also adds a
# `cogmod_ddm5` program (sigmadrift estimated, sigmabias = sigmandt = 0) to
# the tree's benchmarks/gradient_programs.R, since no CI program reaches the
# 5-parameter path.
#
#   Rscript benchmarks/wiener_stan/apply_variant.R --tree <worktree> --variant direct|rec|unrolled|none
#
# `none` adds the cogmod_ddm5 program only (for the base tree). Run from the
# repository root; the tree is meant to be a scratch worktree
# (git worktree add --detach <tmp>/proto HEAD), never this checkout.
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k) { i <- match(paste0("--", k), a); if (is.na(i)) stop("--", k, " missing") else a[i + 1] }
tree <- normalizePath(arg("tree"), winslash = "/")
variant <- arg("variant")
stopifnot(variant == "none" || file.exists(file.path("benchmarks", "wiener_stan", paste0("wfpt_", variant, ".stan"))),
          normalizePath(".", winslash = "/") != tree)

sub1 <- function(x, old, new, fixed = TRUE, perl = FALSE) {
  hits <- gregexpr(old, x, fixed = fixed, perl = perl)[[1]]
  if (sum(hits > 0) != 1) stop("expected exactly one match for: ", substr(old, 1, 60))
  sub(old, new, x, fixed = fixed, perl = perl)
}

if (variant != "none") {
  f <- file.path(tree, "R", "model_ddm.R")
  src <- paste(readLines(f), collapse = "\n")
  stan <- paste(readLines(file.path("benchmarks", "wiener_stan", paste0("wfpt_", variant, ".stan"))),
                collapse = "\n")
  if (grepl("[\"\\\\%]", stan)) stop("the Stan text must hold no quote, backslash or percent sign")
  anchor <- 'formatC(.DDM_TAU0, format = "g", digits = 17, width = 1)), sprintf("\nreal cogmod_ddm_decision_lpdf('
  src <- sub1(src, anchor, paste0('formatC(.DDM_TAU0, format = "g", digits = 17, width = 1)), "\n',
                                  stan, '\n", sprintf("\nreal cogmod_ddm_decision_lpdf('))
  block <- "(?s)  if \\(sw == 0 && sigmandt == 0\\) \\{.*?return wiener_lpdf\\(y \\| boundary, tau0, w, v, sigmadrift\\);\n  \\}"
  src <- sub1(src, block, paste0(
    "  if (sw == 0 && sigmandt == 0) {\n",
    "    // DDM-4 and DDM-5: the Stan-language density of benchmarks/wiener_stan/ (",
    variant, ").\n",
    "    return cogmod_wfpt_lpdf(t | boundary, w, v, sigmadrift);\n",
    "  }"), fixed = FALSE, perl = TRUE)
  writeLines(src, f)
  cat("patched", f, "with", variant, "\n")
}

g <- file.path(tree, "benchmarks", "gradient_programs.R")
gsrc <- paste(readLines(g), collapse = "\n")
if (!grepl("cogmod_ddm5 = list", gsrc, fixed = TRUE)) {
  gsrc <- sub1(gsrc, 'fixed = c(sigmadrift = 0, sigmabias = 0, sigmandt = 0))\n)',
               paste0('fixed = c(sigmadrift = 0, sigmabias = 0, sigmandt = 0)),\n',
                      '  # sigmadrift estimated, the other two variabilities at 0 (the 5-parameter path).\n',
                      '  cogmod_ddm5 = list(family = "cogmod_ddm",\n',
                      '                     fixed = c(sigmabias = 0, sigmandt = 0))\n)'))
  writeLines(gsrc, g)
  cat("added cogmod_ddm5 to", g, "\n")
}
