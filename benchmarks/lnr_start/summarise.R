# Step 2 of the LNR start check (README.md), part 3: the tables, from the
# per-cell CSVs fit_arm.R wrote under results/.
#
#   Rscript benchmarks/lnr_start/summarise.R
#
# Everything is paired by seed: the two arms' chains share their data, program,
# priors and jitter, and differ only in where `Intercept_nuone` starts.
here <- "benchmarks/lnr_start"
R <- file.path(here, "results")
n_warm <- 250L
n_samp <- 250L

read_all <- function(kind) {
  files <- list.files(R, pattern = paste0("__", kind, "\\.csv$"), full.names = TRUE)
  do.call(rbind, lapply(files, function(f) {
    key <- strsplit(sub(paste0("__", kind, "\\.csv$"), "", basename(f)), "__")[[1]]
    cbind(model = key[1], arm = key[2], seed = as.integer(sub("^s", "", key[3])),
          utils::read.csv(f, stringsAsFactors = FALSE))
  }))
}
chains <- read_all("chains")
pars <- read_all("pars")
inits <- read_all("inits")
models <- c("vignette", "vignette_bias", "sim_hier")
models <- models[models %in% chains$model]

# Chains: fit_cell.R keeps a row per chain, complete or not, with its exit
# status and whether cmdstanr counted it (cmdstanr drops a chain whose
# "Elapsed Time" lines it did not parse, exit status regardless) ----------
cat("\nChains complete, of 16 per arm (cmdstanr_dropped: complete CSV that cmdstanr counted as failed)\n")
comp <- do.call(rbind, lapply(split(chains, chains[c("model", "arm")], drop = TRUE),
                              function(x) {
  data.frame(model = x$model[1], arm = x$arm[1], complete = sum(x$complete),
             nonzero_exit = sum(!is.na(x$rc) & x$rc != 0),
             cmdstanr_dropped = sum(x$complete & !x$cmdstanr_ok),
             fits_short = length(unique(x$seed[!x$complete])))
}))
print(comp[order(match(comp$model, models), comp$arm), ], row.names = FALSE)
chains <- chains[chains$complete, ]

# The starts: nuone per arm, and that everything else is the same ----------
cat("\nStarts (Intercept_nuone, link scale; chain 1 of each seed)\n")
st <- inits[inits$par == "Intercept_nuone" & inits$chain == 1, ]
print(stats::reshape(st[, c("model", "seed", "arm", "value")],
                     idvar = c("model", "seed"), timevar = "arm",
                     direction = "wide"), row.names = FALSE)
same <- merge(inits[inits$arm == "A", ], inits[inits$arm == "B", ],
              by = c("model", "seed", "chain", "par"), suffixes = c(".A", ".B"))
same <- same[same$par != "Intercept_nuone", ]
cat(sprintf("Other starts identical between arms: %s (max |A - B| = %.2g)\n",
            all(abs(same$value.A - same$value.B) < 1e-12),
            max(abs(same$value.A - same$value.B))))

# Per fit --------------------------------------------------------------------
fit <- do.call(rbind, lapply(split(chains, chains[c("model", "arm", "seed")],
                                   drop = TRUE), function(x) {
  p <- pars[pars$model == x$model[1] & pars$arm == x$arm[1] & pars$seed == x$seed[1], ]
  data.frame(
    model = x$model[1], arm = x$arm[1], seed = x$seed[1], chains = nrow(x),
    div = sum(x$div), div_warm = sum(x$div_warm), frozen = sum(x$frozen),
    rejected = sum(x$rejected), td10 = sum(x$td10),
    leap_warm = sum(x$leap_warm) / (nrow(x) * n_warm),
    leap_samp = sum(x$leap_samp) / (nrow(x) * n_samp),
    stepsize = stats::median(x$stepsize),
    transient = max(x$transient),
    lp_start = min(x$lp_start),
    wall = max(x$time_total), cpu_samp = sum(x$time_samp),
    ebfmi = min(x$ebfmi),
    max_rhat = max(p$rhat, na.rm = TRUE),
    min_ess = min(p$ess_bulk, na.rm = TRUE),
    med_ess = stats::median(p$ess_bulk, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
fit$min_ess_per_s <- fit$min_ess / fit$wall
fit$med_ess_per_s <- fit$med_ess / fit$wall
fit <- fit[order(fit$model, fit$seed, fit$arm), ]
rownames(fit) <- NULL
utils::write.csv(fit, file.path(R, "fits.csv"), row.names = FALSE)

cat("\nPer arm (over seeds; leapfrog is per draw per chain, transient the worst chain)\n")
arm_tab <- do.call(rbind, lapply(split(fit, fit[c("model", "arm")], drop = TRUE),
                                 function(x) {
  data.frame(
    model = x$model[1], arm = x$arm[1], fits = nrow(x),
    div = sum(x$div), div_warm = sum(x$div_warm), frozen = sum(x$frozen),
    rejected = sum(x$rejected), td10 = sum(x$td10),
    max_rhat = round(max(x$max_rhat), 3),
    leap_warm = round(stats::median(x$leap_warm), 1),
    leap_samp = round(stats::median(x$leap_samp), 1),
    stepsize = signif(stats::median(x$stepsize), 3),
    transient = stats::median(x$transient),
    lp_start = round(stats::median(x$lp_start)),
    min_ess = round(stats::median(x$min_ess)),
    min_ess_per_s = round(stats::median(x$min_ess_per_s), 2),
    wall_s = round(stats::median(x$wall)),
    stringsAsFactors = FALSE
  )
}))
arm_tab <- arm_tab[order(match(arm_tab$model, models), arm_tab$arm), ]
print(arm_tab, row.names = FALSE)

cat("\nPaired ratios B / A by seed, over the seeds where both arms kept all 4 chains (min, median, max over seeds; < 1 is cheaper or worse as marked)\n")
metrics <- c(leap_warm = "cheaper", leap_samp = "cheaper", stepsize = "larger",
             transient = "shorter", wall = "faster", min_ess = "worse",
             min_ess_per_s = "worse", med_ess_per_s = "worse")
full <- fit[fit$chains == 4, ]
ratio_tab <- do.call(rbind, lapply(models, function(m) {
  a <- full[full$model == m & full$arm == "A", ]
  b <- full[full$model == m & full$arm == "B", ]
  s <- intersect(a$seed, b$seed)
  a <- a[match(s, a$seed), ]
  b <- b[match(s, b$seed), ]
  do.call(rbind, lapply(names(metrics), function(k) {
    r <- b[[k]] / a[[k]]
    data.frame(model = m, metric = k, below_1_means = metrics[[k]],
               min = round(min(r), 2), median = round(stats::median(r), 2),
               max = round(max(r), 2), n = length(r), stringsAsFactors = FALSE)
  }))
}))
print(ratio_tab, row.names = FALSE)

# Posterior agreement ----------------------------------------------------------
cat("\nPosterior medians, arm B against arm A, pooled over seeds; z in units of the posterior SD ((q95 - q5) / 3.29)\n")
agree <- do.call(rbind, lapply(models, function(m) {
  p <- pars[pars$model == m, ]
  do.call(rbind, lapply(unique(p$variable), function(v) {
    x <- p[p$variable == v, ]
    a <- x[x$arm == "A", ]
    b <- x[x$arm == "B", ]
    sdv <- stats::median((x$X95. - x$X5.) / 3.29)
    data.frame(model = m, variable = v,
               A = round(stats::median(a$median), 3),
               B = round(stats::median(b$median), 3),
               z = round((stats::median(b$median) - stats::median(a$median)) / sdv, 2),
               max_rhat_A = round(max(a$rhat), 3), max_rhat_B = round(max(b$rhat), 3),
               stringsAsFactors = FALSE)
  }))
}))
print(agree, row.names = FALSE)
utils::write.csv(agree, file.path(R, "agreement.csv"), row.names = FALSE)
