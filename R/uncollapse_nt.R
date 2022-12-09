#' helper to uncollapse nodal types
#'
#' @param nodal_tpyes a List of nodal types
#' @return A List of numeric matrices representing uncollapsed nodal types
#' @keywords internal

uncollapse_nt <- function(nodal_types) {

  x <- nodal_types |>
    lapply(stringr::str_split, "")  |>
    lapply(data.frame) |>
    lapply(t)  |>
    lapply(function(df) apply(df, 2, as.numeric)) |>
    lapply(as.matrix)

  # This is not elegant; to handle cases where a single nodal type exists
  # otherwise it gets wrongly tranposed
  for(j in 1:length(x)) if(length(nodal_types[[j]])==1)  x[[j]] <- t(x[[j]])

  for(j in 1:length(x)){
    # Add row names
    rownames(x[[j]]) <- apply(x[[j]], 1, paste, collapse ="")
    # Add col names
    colnames(x[[j]]) <- CausalQueries:::perm(rep(1, log(ncol(x[[j]]),2))) %>%
      apply(1, paste, collapse = "")
  }

  x
}

