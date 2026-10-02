# The design of the LBA screen: one row per CHAIN, so that every chain is its
# own array task. The question is whether chains of one arm agree, so the unit
# of analysis is the chain, and one chain per task keeps a slow or frozen
# chain from holding a node's worth of others hostage.
#
# Why this screen exists (2026-10-01). The full-data gam_lba fit of the
# Illusion Game project split into two modes (5 chains against 3, 1,300 lp__
# units apart). The cross-check on its siblings found the same in every
# full-data fit that smooths the start-point range: gam_rdm5 (one chain 1,671
# above the other seven) and gam_lnr6 (one chain frozen, one 1,030 below),
# while gam_ddm5, with sigmabias fixed at 0, was clean. Two explanations are on
# the table: the LBA's |v|/s^2 ray, which a condition-varying sigmaone
# multiplies, and the start-point range itself. Re-running the full-data fit
# with the 0.3.4 defaults would mostly show which basin its chains fall into,
# at 4 shards x 2-3 days; this screen asks the question directly, at 200
# participants.
#
#   lba_old   gam_lba, 0.3.3 defaults   does the split exist at this size?
#   lba_new   gam_lba, 0.3.4 defaults   do the defaults change anything?
#   lba_a     sigmaone ~ 1 + (1 | P)    is it the ray?
#   lba_b     sigmabias ~ 1 + (1 | P)   is it the start-point range?
#   rdm5_new  gam_rdm5, 0.3.4 defaults  does the RDM split at this size?
#   rdm5_b    rdm5, sigmabias ~ 1 + P   the same hypothesis, with no ray
#
# "0.3.3 defaults" means the production library's cogmod (d04c7f8, the code
# gam_lba was fitted with); "0.3.4" means this tree. fit_cell.R picks the
# library per arm and checks it got the one it asked for.
#
# Eight chains per arm: with basins taken 5:3, eight chains show a split 98%
# of the time; four would miss it about one time in six.
screen_arms <- function() {
  data.frame(
    arm = c("lba_old", "lba_new", "lba_a", "lba_b", "rdm5_new", "rdm5_b"),
    model = c("gam_lba", "gam_lba", "gam_lba_a", "gam_lba_b", "gam_rdm5",
              "gam_rdm5_b"),
    defaults = c("0.3.3", "0.3.4", "0.3.4", "0.3.4", "0.3.4", "0.3.4"),
    exp = c("lba", "lba", "lba", "lba", "rdm", "rdm"),
    stringsAsFactors = FALSE
  )
}

screen_cells <- function(exp_id) {
  a <- screen_arms()
  if (!exp_id %in% a$exp) stop("unknown experiment '", exp_id, "'")
  a <- a[a$exp == exp_id, ]
  g <- merge(a, data.frame(chain = 1:8), by = NULL)
  g <- g[order(match(g$arm, a$arm), g$chain), ]
  g$id <- seq_len(nrow(g))
  rownames(g) <- NULL
  g[, c("exp", "id", "arm", "model", "defaults", "chain")]
}

# Production settings (analysis/server/fit_model.R): warmup 1000 + 500 draws,
# cold starts. Four threads per chain; every one of these can be overridden
# from the environment for a smoke test (SCREEN_PARTICIPANTS, SCREEN_WARMUP,
# SCREEN_SAMPLES).
screen_layout <- list(participants = 200L, warmup = 1000L, samples = 500L,
                      threads = 4L)
