# Script used to generate data/wagenmakers2008.rda
# Source: https://github.com/DominiqueMakowski/CognitiveModels/blob/main/data/wagenmakers2008.csv

wagenmakers2008 <- read.csv(
    "https://raw.githubusercontent.com/DominiqueMakowski/CognitiveModels/main/data/wagenmakers2008.csv"
)

wagenmakers2008$Participant <- as.integer(wagenmakers2008$Participant)
wagenmakers2008$Condition <- as.character(wagenmakers2008$Condition)
wagenmakers2008$RT <- as.numeric(wagenmakers2008$RT)
# Kept as integer 0/1 rather than logical: `brms` accepts a logical `dec()`
# variable when fitting, but `validate_newdata()` coerces logicals to factors,
# after which `dec()` demands the levels "lower"/"upper" - which makes any
# `newdata` prediction (posterior_predict(), modelbased) impossible.
wagenmakers2008$Error <- as.integer(as.logical(wagenmakers2008$Error))
wagenmakers2008$Frequency <- as.character(wagenmakers2008$Frequency)

save(wagenmakers2008, file = "data/wagenmakers2008.rda", compress = "xz")
