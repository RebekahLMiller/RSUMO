# Define the expected column names for FIMO results
fimo_column_names_cl <-
    c(
        "motif_id",
        "motif_alt_id",
        "sequence_name",
        "start",
        "stop",
        "strand",
        "score",
        "p-value",
        "q-value",
        "matched_sequence"
    )
fimo_column_names_r <-
    c(
        "seqnames",
        "start",
        "end",
        "width",
        "strand",
        "motif_id",
        "motif_alt_id",
        "score",
        "pvalue",
        "qvalue",
        "matched_sequence"
    )

# Define a bunch of dummy variables to get rid of R CMD check notes
motif_id <- motif_alt_id <- values <- ind <- motif_group <- NULL
sequence_name <- seqnames <- matched_sequence <- NULL
score <- TPR <- fg_cumsum <- n_fg <- FPR <- bg_cumsum <- n_bg <- NULL
pvalue <- p_value <- qvalue <- q_value <- NULL
start <- end <- strand <- . <- NULL
