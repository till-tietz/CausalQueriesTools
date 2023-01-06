#' helper to make monotonic non-interacting nodal types for a node with n parents
#' @param n integer specifying number of parents
#' @return a character vector of nodal types
#' @keywords internal
#' @export

simple_nodes <- function(n = 5) {
  N <- 2^n
  s <- 0:(N-1)
  nt <- ((sapply(0:n, function(j) ((s / 2^(j-1)) %% 2 ) >= 1))) |>
    data.frame() |>
    dplyr::mutate(all = 1)

  nt |> apply(2, paste0, collapse = "")
}


#' helper to make model with strictly monotonic relationships and no interactions
#' @param statement string describing causal model
#' @return An object of class \code{causal_model}.
#' @export

make_simple_model <- function(statement){

  b1 <- list(assign(trimws(strsplit(statement,"->|;|<->")[[1]][1]), "01"))
  names(b1) <- trimws(strsplit(statement,"->|;|<->")[[1]][1])
  temp <- CausalQueries::make_model(statement,
                                    nodal_types = b1,
                                    add_causal_types = FALSE)
  n_parents <- temp$dag |>
    dplyr::group_by(children) |>
    dplyr::summarise(n_parents = dplyr::n()) |>
    dplyr::bind_rows(data.frame(children = attr(temp, "root_nodes"), n_parents = 0))
  nodal_types <- lapply(n_parents$n_parents, simple_nodes)
  names(nodal_types) <- n_parents$children
  nodal_types

  m <- CausalQueries::make_model(statement,
                                 nodal_types = nodal_types,
                                 add_causal_types = FALSE)

  m$model_type <- "simple"

  return(m)
}


#' helper to add interactions on selected nodes to a \code{causal_model} made with \link{make_simple_model}
#' @param model \code{causal_model} made with \link{make_simple_model}
#' @param interactions a list of character vectors specifying nodes to interact (the first 1:n-1 nodes per vector specify the interacting nodes, the nth node specifies the interaction target node)
#' @param interaction_operation character vector of logical operations to create interactions with
#' @return An object of class \code{causal_model} with the specified interactions
#' @export

interact_model <- function(model,
                           interactions,
                           interaction_operation = c("&","|")){

  model_type <- model$model_type

  #check model class
  if(!is(model,"causal_model")){
    stop("model must be of class causal_model")
  }

  #check interactions class
  if(!is(interactions,"list")){
    stop("interactions must be a list")
  }

  #check for correct interaction specification i.e. character vectors and  minimum two way interactions on target node
  if(!all(sapply(interactions,is.character))){
    stop("interactions must be specified as character vectors")
  }

  interactions_size <- sapply(interactions,length)

  if(!(all(interactions_size >= 3))){
    stop("each element in interactions should specify a specifc interaction as a vector of node names i.e. two-way interaction of X and Y on Z: \n
         c('X','Y','Z')")
  }

  #check interaction operation specification
  if(!is.character(interaction_operation)){
    stop("interaction_operation should be of type character")
  }

  if(any(!unique(unlist(strsplit(interaction_operation,""))) %in% c("&","|","!"))){
    stop("interaction operations should only be composed of '!','|' or '&'")
  }

  #check if interaction nodes are in model and have parent,child relationship
  error <- ""

  for(i in 1:length(interactions)){

    si <- interactions[[i]]

    target <- si[length(si)]
    inters <- si[1:length(si) - 1]

    nodes_not_in_model <- si[!(si %in% model$nodes)]
    no_relationship <- inters[!(inters %in% model$dag[model$dag$children == target,]$parent)]

    if(length(nodes_not_in_model) != 0){
      nodes_not_in_model <- paste("nodes: ",paste(nodes_not_in_model, collapse = ", "), " not part of model. \n ", sep = "")
    }

    if(length(no_relationship) != 0){
      no_relationship <- paste("no parent relationship between ", paste(no_relationship, collapse = ", "), " and ", target, ".", sep = "")
    }

    if(length(no_relationship) != 0 | length(nodes_not_in_model) != 0){
      error <- paste(error,
                     "\n interactions set ", i, ": c(", paste(interactions[[i]], collapse = ","), "): \n ",
                     nodes_not_in_model,
                     no_relationship," \n \n ",
                     sep = "")
    }

  }

  if(error != ""){
    stop(error)
  }

  #make interaction result vector
  target_nodes <- unique(sapply(interactions, function(i) i[length(i)]))

  res <- vector(mode = "list", length = length(target_nodes))
  names(res) <- target_nodes

  for(i in 1:length(interactions)){

    #get target and interactions
    target <- interactions[[i]][length(interactions[[i]])]
    inters <- interactions[[i]][1:length(interactions[[i]])-1]

    #get monotonicity type for each interaction node on interaction target
    types <- model$nodal_types[[target]][sapply(inters,function(j) which(j == unique(dplyr::filter(model$dag, children == target) |> unlist()))) + 1] |>
      strsplit("") |>
      lapply(function(j) as.logical(as.numeric(j)))

    types <- sapply(interaction_operation, function(j){
      Reduce(j,types) |>
        as.numeric() |>
        as.character() |>
        paste(collapse = "")
    }) |>
      unname()


    res[[target]] <- unique(c(res[[target]],types))
  }

  for(i in 1:length(res)){

    target <- names(res)[i]
    model$nodal_types[[target]] <- unique(c(model$nodal_types[[target]],res[[i]]))

  }

  m <- CausalQueries::make_model(model$statement,
                                 nodal_types = model$nodal_types,
                                 add_causal_types = FALSE)

  if(!is.null(model_type)){
    m$model_type <- paste(model_type,"_interacted", sep = "")
  } else {
    m$model_type <- paste(model_type,"interacted", sep = "")
  }

  m$interactions <- interactions

  return(m)

}

