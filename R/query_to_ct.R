#helper for query deparsing

deparse_query <- function(query, join_by, nodes){

  if (length(query) > 1L){
    stop("Please specify a query of length 1L.")
  }

  if (grepl(".", query, fixed = TRUE)){
    query <- CausalQueries::expand_wildcard(query, join_by = join_by)
  }

  # Global Variables
  i <- 0
  list_names <- ""
  continue <- TRUE

  # strip whitespaces split query into single characters locate opening brackets and reverse oder
  w_query <- gsub(" ", "", query) %>%
    strsplit("") %>%
    unlist()
  bracket_starts <- rev(grep("\\[", w_query))
  bracket_ends <- rev(grep("\\]", w_query))

  if (length(bracket_starts) != length(bracket_ends)) {
    stop("Either '[' or ']' missing.")
  }
  if (length(bracket_starts) == 0) {
    continue = FALSE
  }

  dos <- list()
  vars <- c()
  var_order <- c()

  while (continue) {
    i <- i + 1

    # start at the latest found '[' find the closest subsequent ']' remove brackets and extract
    # expression
    .query <- w_query[(bracket_starts[i]):length(w_query)]
    .bracket_ends <- grep("\\]", .query)[1]
    .query <- .query[1:.bracket_ends]
    brackets <- grepl("\\[|\\]", .query)
    .query <- .query[!brackets]

    # Split expression by ',' 'x = 1, m = 0' into 'x = 1' 'm=0'
    .query <- paste0(.query, collapse = "")
    .query <- unlist(strsplit(.query, ","))

    # Walks through splitted expressions (i.e dos) and evaluates each expression when possible
    if(length(.query) == 0) {
      stop("\nquery does not return any causal types.\nNote that expressions of the form `Y[]==1` are not allowed for mapping queries to causal types.\nSpecify queries as (e.g.) `Y==1` or `Y[X=0] == 1` instead.")
    }

    for (j in 1:length(.query)) {
      do <- paste("list(",.query[j],")", sep = "") %>%
        parse(text = .) %>%
        eval(.,envir = c())

      if (!names(do) %in% nodes){
        stop(paste("Variable", names(do), "is not part of the model."))
      }

      dos <- c(dos,do)
    }

    b <- 1:bracket_starts[i]
    var <- paste0(w_query[b], collapse = "")
    var <- CausalQueries:::st_within(var)
    var <- var[length(var)]
    vars <- c(vars,var)

    # Save result from last iteration and remove corresponding expression w_query
    var_length <- nchar(var)

    .bracket_ends <- bracket_starts[i] + .bracket_ends - 1
    s <- seq(bracket_starts[i] - var_length, .bracket_ends)
    vname <- paste0("var", i)
    var_order <- c(var_order,vname)
    w_query[s[1]] <- vname
    w_query[s[2:length(s)]] <- ""


    # Stop loop there are no [] left
    if (!any(grep("\\[|\\]", w_query))) {
      continue <- FALSE
    }
  }  # end of [] application
  w_query <- w_query[w_query != ""]
  return(list(w_query = w_query,
              dos = rev(dos),
              vars = rev(vars),
              w_query_order = rev(var_order))
         )
}



#map query to causal types
query_to_ct <- function(model,query,join_by = "|", file_name = ""){

  #deparse query
  query_deparsed <- deparse_query(query = query, join_by = join_by, nodes = model$nodes)

  #get number of causal types
  nct <- sapply(model$nodal_types,length) |>
    prod()




}




