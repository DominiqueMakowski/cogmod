# Run benchmarks/gradient_check.R or gradient_cost.R on the DDM-4 (all three
# variabilities fixed at 0 in bf()), which the CI programs do not build:
# gradient_programs.R estimates every DDM dpar, so CI only ever reaches the
# 7-parameter call site. The harness is left untouched; the script is
# evaluated with .GP_FIXED extended right after gradient_programs.R loads.
#
#   Rscript benchmarks/ddm_router/ddm4_harness.R benchmarks/gradient_check.R --families cogmod_ddm --out DIR
#   Rscript benchmarks/ddm_router/ddm4_harness.R benchmarks/gradient_cost.R emit --pkg <base> --out DIR --families cogmod_ddm
args <- commandArgs(trailingOnly = TRUE)
code <- readLines(args[1])
hook <- 'source("benchmarks/gradient_programs.R")'
i <- which(startsWith(trimws(code), hook))
stopifnot(length(i) >= 1)
code[i[1]] <- paste0(code[i[1]], '; .GP_FIXED$cogmod_ddm <- c(sigmadrift = 0, sigmabias = 0, sigmandt = 0)')
# commandArgs(trailingOnly = TRUE) inside the script must see only its own
# options; gradient_cost.R wraps base::commandArgs, so point it at this one.
code <- sub("orig <- base::commandArgs", "orig <- commandArgs", code, fixed = TRUE)
rest <- args[-1]
commandArgs <- function(trailingOnly = FALSE) if (trailingOnly) rest else base::commandArgs(FALSE)
eval(parse(text = code), envir = globalenv())
