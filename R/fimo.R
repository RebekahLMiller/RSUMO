# #' Run FIMO
# #'
# #' Run FIMO on foreground and background sequences with given motifs
# #'
# #' @param motif_library meme file
# #' @param fg_sequences fasta file
# #' @param bg_sequences fasta file or NA
# #'
# #' @return a data frame of fimo results
# #' @export
# #'
# #' @examples
# #' print("")
# run_fimo <- function(motif_library, fg_sequences, bg_sequences = NA) {
#     # Get a list of all the motif names
#     motif_names <-
#         universalmotif::read_meme(motif_library) %>%
#
#         # Convert each motif to a data frame
#         lapply(universalmotif::as.data.frame) %>%
#
#         # Combine all the data frames
#         dplyr::bind_rows() %>%
#
#         # Keep just the motif name and altname columns
#         dplyr::select(motif_id = name, motif_alt_id = altname)
#
#     # Get a list of all the foreground sequence names
#     fg_seq_names <-
#         Biostrings::readDNAStringSet(fg_sequences) %>%
#         names()
#
#     # Get a list of all the background sequence names
#     if (!is.na(bg_sequences)) {
#         bg_seq_names <-
#             Biostrings::readDNAStringSet(bg_sequences) %>%
#             names()
#     }
#
#     # Run FIMO on the foreground sequences
#     fg_fimo_results <-
#         memes::runFimo(
#             sequences = fg_sequences,
#             motifs = motif_library,
#             parse_genomic_coord = FALSE,
#             skip_matched_sequence = TRUE,
#             thresh = 0.01
#         ) %>%
#
#         # Convert to a data frame
#         as.data.frame() %>%
#
#         # Convert the sequence names and motif names into factors
#         dplyr::mutate(
#             seqnames = factor(seqnames, levels = fg_seq_names),
#             motif_id = factor(motif_id, levels = motif_names$motif_id)
#         ) %>%
#
#         # Make sure all combinations of sequences and motifs are present
#         tidyr::complete(
#             seqnames,
#             motif_id,
#             fill = list(score = 0)
#         )
#
#     # Run FIMO on the background sequences
#     if (!is.na(bg_sequences)) {
#         bg_fimo_results <-
#             memes::runFimo(
#                 sequences = bg_sequences,
#                 motifs = motif_library,
#                 parse_genomic_coord = FALSE,
#                 skip_matched_sequence = TRUE,
#                 thresh = 0.01
#             ) %>%
#
#             # Convert to a data frame
#             as.data.frame() %>%
#
#             # Convert the sequence names and motif names into factors
#             dplyr::mutate(
#                 seqnames = factor(seqnames, levels = bg_seq_names),
#                 motif_id = factor(motif_id, levels = motif_names$motif_id),
#             ) %>%
#
#             # Make sure all combinations of sequences and motifs are present
#             tidyr::complete(
#                 seqnames,
#                 motif_id,
#                 fill = list(score = 0)
#             )
#
#         # Combine the results from the foreground and background
#         all_fimo_results <-
#             dplyr::bind_rows(
#                 "fg" = fg_fimo_results,
#                 "bg" = bg_fimo_results,
#                 .id = "source"
#             )
#     } else {
#         # If there are no bg sequences, add a source column to the fg table
#         all_fimo_results <-
#             dplyr::mutate(fg_fimo_results, source = "fg", .before = seqnames)
#     }
#
#     # Fill in missing motif alt names
#     all_fimo_results <-
#          all_fimo_results %>%
#
#         # Get missing motif alt names
#         dplyr::left_join(
#             motif_names,
#             by = "motif_id"
#         ) %>%
#
#         # Combine the two motif alt name columns into one
#         dplyr::mutate(
#             motif_alt_id = dplyr::coalesce(motif_alt_id.x, motif_alt_id.y),
#             .keep = "unused",
#             .after = motif_id
#         )
#
#     # Return the table of all FIMO results
#     return(all_fimo_results)
# }

#' Preprocess FIMO results
#'
#' Prepare FIMO results for downstream analyses by filling in any missing
#' combinations of sequences and motifs and reformatting the table if necessary.
#'
#' @param fimo_results `data.frame`. The table of FIMO results to process.
#' @param sequences `character` or `NA`. The list of all sequences that should
#'   be in the final table. If `NA`, all of the sequences in the input
#'   `fimo_results` table that have at least one motif that passes the p-value
#'   threshold are used. (Default: NA)
#' @param motifs `character` or `NA`. The list of all motifs that should be in
#'   the final table. If `NA`, all of the motifs in the input `fimo_results`
#'   table that have at least one occurrence that passes the p-value threshold
#'   are used. (Default: NA)
#' @param max_pvalue `numeric(1)`. The maximum p-value to consider. Set to `1`
#'   to skip filtering by p-value. (Default: 1e-4)
#'
#' @return A `data.frame` with columns `motif_id`, `motif_alt_id`,
#'   `sequence_name`, `start`, `stop`, `strand`, `score`, `p_value`, `q_value`,
#'   and `matched_sequence`.
#' @export
#'
#' @examples
#' print("")
preprocess_fimo_results <-
    function(
        fimo_results,
        sequences = NA,
        motifs = NA,
        max_pvalue = 1e-4
    ) {
    # Make sure the FIMO results are a data frame
    fimo_results <- as.data.frame(fimo_results)

    # Make sure the column names are one of the two expected possibilities
    if (identical(colnames(fimo_results), fimo_column_names_cl)) {
        # Replace "-"s with "_"s in command line output
        fimo_results <-
            dplyr::rename(
                fimo_results,
                p_value = "p-value",
                q_value = "q-value"
            )
    } else if (identical(colnames(fimo_results), fimo_column_names_r)) {
        # Make the R wrapper output consistent with the command line output
        fimo_results <-
            dplyr::select(
                fimo_results,
                motif_id:motif_alt_id,
                sequence_name = seqnames,
                start,
                stop = end,
                strand,
                score,
                p_value = pvalue,
                q_value = qvalue,
                matched_sequence
            )
    } else {
        stop(
            "Unexpected column names: where did you get these FIMO results? :(",
            call. = FALSE
        )
    }

    # Keep only FIMO results that pass the p-value threshold
    fimo_results <-
        dplyr::filter(fimo_results, p_value <= max_pvalue)

    # If no sequence names are given, use all the sequences in the FIMO results
    if (is.na(sequences)) {
        sequences <- unique(fimo_results$sequence_name)
    } else {
        fimo_results <-
            dplyr::filter(fimo_results, sequence_name %in% sequences)
    }

    # If no motif names are given, use all the motifs in the FIMO results
    if (is.na(motifs)) {
        motifs <- unique(fimo_results$motif_id)
    } else {
        fimo_results <-
            dplyr::filter(fimo_results, motif_id %in% motifs)
    }

    # Fill in any missing combinations of sequences and motifs
    fimo_results <-
        fimo_results %>%

        # Make the sequences and motifs factors and the q-values numeric
        dplyr::mutate(
            motif_id = factor(motif_id, levels = motifs),
            sequence_name = factor(sequence_name, levels = sequences),
            q_value = as.numeric(q_value)
        ) %>%

        # Make sure all combinations of sequences and motifs are present
        tidyr::complete(
            tidyr::nesting(motif_id, motif_alt_id),
            sequence_name,
            fill =
                list(
                    score = 0,
                    p_value = 1,
                    q_value = 1
                )
        ) %>%

        # Convert the sequence names and motifs names back into not factors
        dplyr::mutate(
            motif_id = as.character(motif_id),
            sequence_name = as.character(sequence_name)
        )
}

