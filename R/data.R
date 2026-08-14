#' Lexical Decision Data from Wagenmakers et al. (2008)
#'
#' Response times and accuracy from a lexical decision task, Experiment 1 in
#' \href{https://doi.org/10.1016/j.jml.2007.04.006}{Wagenmakers et al. (2008)}
#' (also reanalyzed in
#' \href{https://doi.org/10.3389/fpsyg.2012.00292}{Heathcote & Love, 2012}).
#' Six participants judged whether letter strings were valid English words
#' under instructions that emphasized either speed or accuracy.
#'
#' @format A data frame with 10,369 rows and 5 variables:
#' \describe{
#'   \item{Participant}{Participant identifier (integer, 1-6).}
#'   \item{Condition}{Instruction condition, `"Speed"` or `"Accuracy"`.}
#'   \item{RT}{Response time, in seconds.}
#'   \item{Error}{Logical, whether the response was an error.}
#'   \item{Frequency}{Word frequency category of the presented stimulus.}
#' }
#'
#' @source
#' <https://github.com/DominiqueMakowski/CognitiveModels/blob/main/data/wagenmakers2008.csv>
#'
#' @examples
#' data(wagenmakers2008)
#' head(wagenmakers2008)
"wagenmakers2008"


#' Simulated Data Where Linear Models Fail
#'
#' A simulated repeated-measures experiment in which the two conditions have,
#' **by construction, the same mean reaction time** while differing radically in
#' every other respect. Condition `A` has a short non-decision time (50 ms) and
#' a wide, heavily right-skewed distribution; condition `B` has a long
#' non-decision time (450 ms) and a narrow, nearly symmetric one. Twenty
#' participants each contribute 25 trials per condition.
#'
#' The dataset exists to demonstrate the limits of the summary-statistics
#' approach: a linear model - including a correctly specified linear *mixed*
#' model with a random intercept per participant - finds no effect of
#' `Condition`, because a difference in shift, in spread and in tail weight is
#' invisible to a comparison of means. It is used in the *RT models* vignette
#' and in the `cogmod` paper.
#'
#' Participants differ in their overall speed through an *additive* offset
#' (SD = 30 ms) applied identically to both conditions, so that each
#' participant's true condition effect is exactly zero and a random intercept is
#' the correctly specified model for the between-participant variation. The
#' latent offsets are not included in the data.
#'
#' @format A data frame with 1,000 rows and 3 variables:
#' \describe{
#'   \item{Participant}{Participant identifier (factor, `S01`-`S20`).}
#'   \item{Condition}{Experimental condition, `"A"` or `"B"`.}
#'   \item{RT}{Reaction time, in seconds.}
#' }
#'
#' @source Simulated; see `data-raw/badlm.R` for the generating code.
#'
#' @examples
#' data(badlm)
#'
#' # The two conditions have the same mean...
#' tapply(badlm$RT, badlm$Condition, mean)
#'
#' # ...but nothing else in common
#' tapply(badlm$RT, badlm$Condition, sd)
#'
#' \donttest{
#' # A mixed model finds nothing
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   summary(lme4::lmer(RT ~ Condition + (1 | Participant), data = badlm))
#' }
#' }
"badlm"
