# Downloads the ds000030 task-switching event files used by example_cnp.R.
# Split out from example_cnp.R so the (slow, network-bound) download can be
# run once and independently of the model fit.
#
#   Rscript download_cnp.R
#
# Run from the `paper/` directory. Writes cnp_data/participants.tsv and
# cnp_data/events/<sub>.tsv. Re-running skips files already on disk.

BASE <- "https://s3.amazonaws.com/openneuro.org/ds000030"
DIR <- "cnp_data"
dir.create(file.path(DIR, "events"), recursive = TRUE, showWarnings = FALSE)

pfile <- file.path(DIR, "participants.tsv")
if (!file.exists(pfile)) download.file(file.path(BASE, "participants.tsv"), pfile, quiet = TRUE)
participants <- read.delim(pfile)

subs <- participants$participant_id[
  participants$diagnosis %in% c("CONTROL", "ADHD") & participants$taskswitch == "1"
]
cat(sprintf("%d candidate participants\n", length(subs)))

ok <- 0
for (i in seq_along(subs)) {
  s <- subs[i]
  out <- file.path(DIR, "events", paste0(s, ".tsv"))
  if (file.exists(out) && file.size(out) > 500) { ok <- ok + 1; next }
  url <- sprintf("%s/%s/func/%s_task-taskswitch_events.tsv", BASE, s, s)
  res <- try(download.file(url, out, quiet = TRUE), silent = TRUE)
  if (!inherits(res, "try-error") && file.exists(out) && file.size(out) > 500) {
    ok <- ok + 1
  } else if (file.exists(out)) {
    unlink(out)
  }
  if (i %% 20 == 0) cat(sprintf("  %d/%d checked, %d files\n", i, length(subs), ok))
}
cat(sprintf("done: %d event files\n", ok))
