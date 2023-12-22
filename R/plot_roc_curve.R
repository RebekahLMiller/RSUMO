#' Calculate AUC
#'
#' Calculate the AUC of the ROC curve
#'
#' @param ranked_scores a data frame of scored sequences with FPR and TPR
#'   columns
#' @param alternative `character(1)`. One of "two.sided, "less", or "greater".
#'   The alternative hypothesis to test.
#'
#' @return a `data.frame` with columns `motif_id`, `auc`, `w_statistic`,
#'   `p_value`, and `alternative.`
#' @export
#'
#' @examples
#' print("")
calculate_auc <-
    function(
        ranked_scores,
        alternative = c("two.sided", "less", "greater")
    ) {
    # Make sure the selected alternative hypothesis is a valid option
    alternative <- match.arg(alternative)

    # Get a list of all the motifs to calculate the AUC for
    motif_ids <- unique(ranked_scores$motif_id)

    # Calculate the AUC for each motif
    aucs <-
        lapply(motif_ids, function(motif) {
            # Filter the table to just rows for this motif
            ranked_scores_motif <-
                dplyr::filter(ranked_scores, motif_id == motif)

            # Pull out the foreground scores
            fg_scores <-
                dplyr::filter(ranked_scores_motif, source == "fg") %>%
                dplyr::pull(score)

            # Pull out the background scores
            bg_scores <-
                dplyr::filter(ranked_scores_motif, source == "bg") %>%
                dplyr::pull(score)

            # Perform the Wilcoxon rank sum test
            test_result <-
                stats::wilcox.test(
                    fg_scores,
                    bg_scores,
                    alternative = alternative
                )

            # Calculate the AUC
            auc <-
                test_result$statistic / (length(fg_scores) * length(bg_scores))

            # Make a data frame with all the useful values
            data.frame(
                "motif_id" = motif,
                "auc" = auc,
                "w_statistic" = test_result$statistic,
                "p_value" = test_result$p.value,
                "alternative" = test_result$alternative,
                row.names = NULL
            )
        }) %>%

        # Combine all the data frames into one
        dplyr::bind_rows()

    # Return the data frame of AUCs
    return(aucs)
}

#' Plot ROC curves
#'
#' Plot a ROC curve
#'
#' @param ranked_scores a data frame of scored sequences with FPR and TPR
#'   columns
#'
#' @return a ggplot object with the ROC curve
#' @export
#'
#' @examples
#' print("")
plot_roc_curve <- function(ranked_scores) {
    # Handle empty data frames
    if (nrow(ranked_scores) == 0) {
        return(ggplot2::ggplot() + ggplot2::geom_blank())
    }

    # Calculate the AUC
    auc <- calculate_auc(ranked_scores)

    # Plot the ROC curve
    roc_curve_plot <-
        ggplot2::ggplot(ranked_scores, ggplot2::aes(x = FPR, y = TPR)) +

        # Plot the curve
        ggplot2::geom_line(color = "red4") +

        # Add a diagonal line
        ggplot2::geom_abline(
            intercept = 0,
            slope = 1,
            color = "grey",
            linetype = "dashed"
        ) +

        # Add the AUC and p-value
        ggplot2::annotate(
            "label",
            x = 0.85,
            y = 0.15,
            label = paste(
                paste("AUC:", signif(auc$auc, 3)),
                paste("p-value:", signif(auc$p_value, 3)),
                paste("alternative:", auc$alternative),
                sep = "\n"
            )
        ) +

        # Change the overall theme
        ggplot2::theme_bw()

    # Return the plot
    return(roc_curve_plot)
}

