#' extract variables and associated do-operations from query
#' and re-write query in terms of realized outcomes post do-operation
#'
#' @param query String specifying query
#' @param join_by Logical operator. Used to connect causal statements: \emph{AND} ('&') or \emph{OR} ('|').
#' @param nodes Character String of nodes in model
#' @return A List containing, the reconstructed query, a list of do-operations, a vector of node names attached to the dos, a vector of newly constructed variables
#' @keywords internal
#' @export


deparse_query <- function(query, join_by, nodes){

  if (length(query) > 1L){
    stop("Please specify a query of length 1L.")
  }

  if (grepl(".", query, fixed = TRUE)){
    query <- CausalQueries::expand_wildcard(query, join_by = join_by)
  }

  #split query into individual characters
  w_query <- gsub(" ", "", query) |>
    strsplit("") |>
    unlist()

  #detect double character operators (==,>=,>=) and re-join
  double_operator <- grep(">|<|=",w_query)
  double_operator_id <- which(diff(double_operator) == 1)

  if(length(double_operator_id) > 0){
    for(i in double_operator_id){
      w_query[double_operator[i]] <- paste(w_query[double_operator[i]],w_query[double_operator[i]+1], sep = "")
    }
    w_query <- w_query[-(double_operator[double_operator_id]+1)]
  }

  #re-join variable names (find non variable name symbols, pad with white space,
  #collapse vector into string, split at white space)
  sym <- paste(c("\\[","\\]","=","==",">=",">","<","<=",
                 "\\&","\\|",",","\\+","\\-"),collapse = "|")
  sym_id <- grep(sym,w_query)

  w_query[sym_id] <- paste0(" ",w_query[sym_id]," ")
  w_query <- paste(w_query,collapse = "")|>
    strsplit(" ")|>
    unlist()
  w_query <- w_query[w_query != ""]

  #split query into component sets
  any(grepl("\\&", w_query)) {

  } else {

  }


  set_operation <- which(head(w_query, -1) == "]" & tail(w_query, -1) %in% c("==",">=","<=",">","<","+","-")) + 1
  w_query_1 <- w_query[1:(set_operation-1)]
  w_query_2 <- w_query[(set_operation+1):length(w_query)]

}



deparse_query <- function(query, join_by, nodes){

  if (length(query) > 1L){
    stop("Please specify a query of length 1L.")
  }

  if (grepl(".", query, fixed = TRUE)){
    query <- CausalQueries::expand_wildcard(query, join_by = join_by)
  }

  #split query into individual characters
  w_query <- gsub(" ", "", query) |>
    strsplit("") |>
    unlist()

  #detect double character operators (==,>=,>=) and re-join
  double_operator <- grep(">|<|=",w_query)
  double_operator_id <- which(diff(double_operator) == 1)

  if(length(double_operator_id) > 0){
    for(i in double_operator_id){
      w_query[double_operator[i]] <- paste(w_query[double_operator[i]],w_query[double_operator[i]+1], sep = "")
    }
    w_query <- w_query[-(double_operator[double_operator_id]+1)]
  }

  #re-join variable names (find non variable name symbols, pad with white space,
  #collapse vector into string, split at white space)
  sym <- paste(c("\\[","\\]","=","==",">=",">","<","<=",
                 "\\&","\\|"),collapse = "|")
  sym_id <- grep(sym,w_query)

  w_query[sym_id] <- paste0(" ",w_query[sym_id]," ")
  w_query <- paste(w_query,collapse = "")|>
    strsplit(" ")|>
    unlist()
  w_query <- w_query[w_query != ""]

  #find positions of variables and brackets
  node_pos <- grep(paste(nodes, collapse = "|"),w_query)
  bracket_starts <- rev(grep("\\[", w_query))
  bracket_ends <- rev(grep("\\]", w_query))

  if (length(bracket_starts) != length(bracket_ends)) {
    stop("Either '[' or ']' missing.")
  }

  #drop variables within brackets from bracket positions
  drop <- sapply(node_pos, function(i){
    sapply(1:length(bracket_starts), function(j){
      bracket_starts[j] < i && i < bracket_ends[j]
    })|>
      any()
  })

  node_pos <- rev(node_pos[!drop])
  dos <- vector(mode = "list", length = length(node_pos))
  vars <- c()
  var_order <- c()

  for(i in 1:length(node_pos)){

    if(w_query[node_pos[i] + 1] != "["){
      do <- list()
      do_j <- paste("list(",w_query[node_pos[i]]," = -1 )", sep = "") |>
        parse(file = "", n = NULL) |>
        eval(envir = c())
      do <- c(do,do_j)
      dos[[i]] <- do
      vars <- c(vars,w_query[node_pos[i]])
      vname <- paste("var",i,sep="")
      var_order <- c(var_order,vname)
      w_query[node_pos[i]] <- vname
    } else {
      open_bracket <- node_pos[i] + 2
      close_bracket <- bracket_ends[bracket_ends > open_bracket]
      close_bracket <- close_bracket[which.max(open_bracket - close_bracket)] - 1
      sub_query <- w_query[open_bracket:close_bracket]

      # Split expression by ',' 'x = 1, m = 0' into 'x = 1' 'm=0'
      sub_query <- paste0(sub_query, collapse = "")
      sub_query <- unlist(strsplit(sub_query, ","))

      # Walks through splitted expressions (i.e dos) and evaluates each expression when possible
      if(length(sub_query) == 0) {
        stop("\nquery does not return any causal types.\nNote that expressions of the form `Y[]==1` are not allowed for mapping queries to causal types.\nSpecify queries as (e.g.) `Y==1` or `Y[X=0] == 1` instead.")
      }

      do <- list()
      for (j in 1:length(sub_query)) {
        do_j <- paste("list(",sub_query[j],")", sep = "") |>
          parse(file = "", n = NULL) |>
          eval(envir = c())

        if (!names(do_j) %in% nodes){
          stop(paste("Variable", names(do_j), "is not part of the model."))
        }
        do <- c(do,do_j)
      }
      dos[[i]] <- do
      vars <- c(vars,w_query[node_pos[i]])
      vname <- paste("var",i,sep="")
      var_order <- c(var_order,vname)
      w_query[node_pos[i]] <- vname
      w_query <- w_query[-((open_bracket - 1):(close_bracket + 1))]
    }
  }

  #turn w_query into list of pairwise operations
  ind <- c(1,grep("\\&",w_query),length(w_query))
  w_query <- lapply(1:(length(ind)-1), function(i){
    string <- w_query[ind[i]:ind[i+1]]
    string <- string[string != "&"]
  })

  dos[sapply(dos, is.null)] <- NULL
  return(list(w_query = w_query,
              dos = rev(dos),
              vars = rev(vars),
              w_query_order = rev(var_order))
  )
}







