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
