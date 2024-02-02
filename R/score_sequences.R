#' Score sequences
#'
#' Calculate a single score for each sequence for each individual motif.
#'
#' @param fimo_results `data.frame`. The table of FIMO results to use. See
#'   [preprocess_fimo_results()] for expected columns.
#' @param method `character(1)`. One of "sum_pvalue", "max_pvalue", "sum_score",
#'   or "max_score". The scoring method to use. (Default: "sum_pvalue")
#'
#' @return A `data.frame` with columns `sequence_name`, `motif_id`, and `score`.
#' @export
#'
#' @examples
#' print("")
score_sequences <-
    function(
        fimo_results,
        method = c("sum_pvalue", "max_pvalue", "sum_score", "max_score")
    ) {
    # Make sure the selected scoring method is a valid option
    method <- match.arg(method)

    # Calculate the appropriate score for each sequence and each motif
    scored_sequences <-
        data.table::setDT(fimo_results)[
            ,
            .(
                score = dplyr::case_when(
                    method == "sum_pvalue" ~ sum(-log10(p_value)),
                    method == "max_pvalue" ~ max(-log10(p_value)),
                    method == "sum_score" ~ sum(score[score > 0]),
                    method == "max_score" ~ max(score)
                )
            ),
            by = .(sequence_name, motif_id)
        ]

    # Return the table of sequence scores
    return(scored_sequences)
}

#' Combine motif scores
#'
#' Calculate a single score for each sequence based on multiple motifs.
#'
#' @param scored_sequences `data.frame`. The table of sequence scores to
#'   combine.
#' @param motif_groups `list`. A named list of motif groupings to score.
#' @param method `character(1)`. One of "sum", "max", or "mean". The scoring
#'   method to use. (Default: "sum")
#'
#' @return A `data.frame` with columns `sequence_name`, `motif_id`, and `score`,
#'   where the IDs in the  `motif_id` column are the names of the motif
#'   groupings given in `motif_groups`.
#' @export
#'
#' @examples
#' print("")
combine_scores <-
    function(
        scored_sequences,
        motif_groups,
        method = c("sum", "max", "mean")
    ) {
    # Make sure the selected scoring method is a valid option
    method <- match.arg(method)

    # Make sure each group has a unique name
    if (length(unique(names(motif_groups))) != length(motif_groups)) {
        stop(
            "Motif groups do not have unique names",
            call. = FALSE
        )
    }

    # Make a data frame of the motif groupings
    motif_groups <-
        utils::stack(motif_groups) %>%

        # Rename the columns
        dplyr::rename(
            motif_id = values,
            motif_group = ind
        )

    # Add a column with the motif groupings
    scored_sequences <-
        dplyr::inner_join(
            scored_sequences,
            motif_groups,
            by = "motif_id"
        )

    # Calculate the appropriate score for each sequence and each motif
    scored_sequences_grouped <-
        data.table::setDT(scored_sequences)[
            ,
            .(
                score = dplyr::case_when(
                    method == "sum" ~ sum(score),
                    method == "max" ~ max(score),
                    method == "mean" ~ mean(score)
                )
            ),
            by = .(sequence_name, motif_group)
        ]

    # Rename the motif_group column
    scored_sequences_grouped <-
        dplyr::rename(scored_sequences_grouped, motif_id = motif_group)

    # Return the table of sequence scores
    return(scored_sequences_grouped)
}

#' Calculate TPR and FPR
#'
#' The title is good enough for now
#'
#' @param scored_sequences_fg `data.frame`. Sequences scored by motif
#'   occurrences.
#' @param scored_sequences_bg `data.frame`. Sequences scored by motif
#'   occurrences.
#'
#' @return A `data.frame` of scored sequences with FPR and TPR columns added
#' @export
#'
#' @examples
#' print("")
rank_scores <- function(scored_sequences_fg, scored_sequences_bg) {
    # Rank the scores for each motif and calculate the TPR and FPR
    ranked_scores <-
        dplyr::bind_rows(
            "fg" = scored_sequences_fg,
            "bg" = scored_sequences_bg,
            .id = "source"
        ) %>%

        # Shuffle so that after sorting by scores, it's not always fg then bg
        dplyr::slice_sample(., n = nrow(.), replace = FALSE) %>%

        # Group by motif
        dplyr::group_by(motif_id) %>%

        # Sort by score within each group
        dplyr::arrange(dplyr::desc(score), .by_group = TRUE) %>%

        # Add cumulative counts of foreground/background sequences and TPR/FPR
        dplyr::mutate(
            fg_cumsum = cumsum(source == "fg"),
            bg_cumsum = cumsum(source == "bg"),
            n_fg = sum(source == "fg"),
            n_bg = sum(source == "bg"),
            TPR = fg_cumsum / n_fg,
            FPR = bg_cumsum / n_bg
        ) %>%

        # Remove the counts of foreground/background sequences
        dplyr::select(-n_fg, -n_bg) %>%

        # Remove grouping
        dplyr::ungroup()

    # Return the ranked scores
    return(ranked_scores)
}

#' Get genes with motifs
#'
#' Get a list of genes that contain one or more given motifs
#'
#' @param scored_sequences a data frame of FIMO results
#' @param motif_ids the IDs of the motifs to look for
#' @param min_score the minimum score to count a motif occurrence
#'
#' @return a data frame of sequences that have the given motifs
#' @export
#'
#' @examples
#' print("")
get_sequences_with_motifs <-
    function(
        scored_sequences,
        motif_ids,
        min_score = 1
    ) {
        # Get a list of names of sequences that have a given motif
        sequences <-
            scored_sequences %>%

            # Keep only rows that have a given motif with a high enough score
            dplyr::filter(
                motif_id %in% motif_ids,
                score >= min_score
            ) %>%

            # Pull out the sequence names
            dplyr::pull(sequence_name) %>%

            # Remove duplicates
            unique()

        # Return the list of sequences that have a given motif
        return(sequences)
    }

