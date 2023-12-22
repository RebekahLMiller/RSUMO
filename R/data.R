#' Example table of sequence scores
#'
#' A table containing 1,348 promoter sequences scored for 10 motifs.
#'
#' @format A data frame with 13,480 rows and 3 columns:
#' \describe{
#'   \item{sequence_name}{The name of the promoter sequence.}
#'   \item{motif_id}{The name of the motif.}
#'   \item{score}{The score for this motif in this promoter sequence.}
#' }
"example_scored_sequences"

#' Example table of ranked sequence scores
#'
#' A table containing the rankings of 590 foreground promoter sequences and 758
#' background promoter sequences scored for 10 motifs.
#'
#' @format A data frame with 13,480 rows and 8 columns:
#' \describe{
#'   \item{source}{Whether this sequence is part of the foreground ("fg") set
#'    or the background ("bg") set.}
#'   \item{sequence_name}{The name of the promoter sequence.}
#'   \item{motif_id}{The name of the motif.}
#'   \item{score}{The score for this motif in this promoter sequence.}
#'   \item{fg_cumsum}{The number of foreground sequences ranked as well as or
#'   better than this sequence for this motif.}
#'   \item{bg_cumsum}{The number of background sequences ranked as well as or
#'   better than this sequence for this motif.}
#'   \item{TPR}{The fraction of foreground sequences ranked as well as or
#'   better than this sequence for this motif.}
#'   \item{FPR}{The fraction of background sequences ranked as well as or
#'   better than this sequence for this motif.}
#' }
"example_ranked_scores"

