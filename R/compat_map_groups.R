# Compatibility helper for legacy `.map_groups_dfr()` calls --------------------

if (!exists(".map_groups_dfr", mode = "function", inherits = TRUE)) {
  .map_groups_dfr <- function(.data, .f, ..., .drop_keys = TRUE, .template = NULL) {
    if (!dplyr::is_grouped_df(.data)) {
      stop("`.map_groups_dfr()` requires a grouped data frame.")
    }

    fn <- rlang::as_function(.f)

    pieces <- dplyr::group_map(
      .data,
      function(.x, .key) {
        out <- fn(.x, .key, ...)

        if (is.null(out)) {
          out <- tibble::tibble()
        } else {
          out <- tibble::as_tibble(out)
        }

        if (!is.null(.template) && !nrow(out)) {
          out <- .template[0, , drop = FALSE]
        }

        if (.drop_keys) {
          keep_cols <- setdiff(names(out), names(.key))
          out <- out[, keep_cols, drop = FALSE]

          if (!nrow(out)) {
            return(out)
          }

          key_df <- .key[rep(1, nrow(out)), , drop = FALSE]
          return(dplyr::bind_cols(key_df, out))
        }

        if (!nrow(out)) {
          return(out)
        }

        missing_keys <- setdiff(names(.key), names(out))
        if (length(missing_keys)) {
          key_df <- .key[rep(1, nrow(out)), missing_keys, drop = FALSE]
          out <- dplyr::bind_cols(key_df, out)
        }

        out
      },
      .keep = TRUE
    )

    if (!length(pieces)) {
      if (!is.null(.template)) {
        return(.template[0, , drop = FALSE])
      }
      return(tibble::tibble())
    }

    res <- dplyr::bind_rows(pieces)
    if (!nrow(res) && !is.null(.template)) {
      return(.template[0, , drop = FALSE])
    }

    res
  }
}

