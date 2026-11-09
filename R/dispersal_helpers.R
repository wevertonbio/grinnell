# dispersers according to values from suitability layer
disperser_rules <- function(layer, max_dispersers) {
  maxsuit <- terra::minmax(layer)[2]
  if (!is.finite(maxsuit)) maxsuit <- 0
  Nd_list <- seq(1, max_dispersers)
  incNd <- maxsuit / max_dispersers
  S_list <- seq(incNd, maxsuit, incNd)

  return(list(Nd_list = Nd_list, S_list = S_list))
}

# defines number of dispersers according to suitability and rules
nd_sval <- function(S_values, disperser_rules) {
  S_list <- disperser_rules$S_list
  Nd_list <- disperser_rules$Nd_list
  idx <- findInterval(S_values, S_list, rightmost.closed = TRUE) + 1
  idx[idx > length(Nd_list)] <- length(Nd_list)
  return(Nd_list[idx])
}


# finds locations in the matrix according to coordinates and matrix features
set_loc <- function(coord, NW_vertex, cell_size) {
  loc_val <- cbind(row = NW_vertex[1] - coord[, 2],
                   column = coord[, 1] - NW_vertex[2])

  return(round(loc_val / cell_size))
}


# creates a matrix where all cells with records are set as 1 and the others as 0
set_pop <- function(data, long, lat, NW_vertex, layer_dim, cell_size,
                    proportion = 1, rule = "random", set_seed = 1,
                    suit_values = NULL) {
  if (proportion == 1) {
    samp <- data[, c(long, lat)]
  } else {
    nr <- nrow(data)
    n <- round(nr * proportion)

    set.seed(set_seed)
    if (rule == "random") {
      samp <- data[sample(nr, n), c(long, lat)]
    } else {
      if (is.null(suit_values)) {
        stop("Argument 'suit_values' must be defined when rule = 'suitability'")
      }
      suit_values <- as.numeric(suit_values)
      suit_values[is.na(suit_values)] <- 0
      if (length(suit_values) != nr) {
        stop("'suit_values' must have one value per row in 'data'")
      }
      if (sum(suit_values) <= 0) {
        stop("All suitability values are <= 0; cannot sample by suitability")
      }
      samp <- data[sample(nr, n, prob = suit_values), c(long, lat)]
    }
  }

  cells <- set_loc(samp, NW_vertex, cell_size)
  C <- matrix(data = 0, nrow = layer_dim[1], ncol = layer_dim[2])

  # linear indices are faster than assigning row/column pairs one by one
  lin <- (cells[, 2] - 1L) * layer_dim[1] + cells[, 1]
  valid <- !is.na(lin) & lin >= 1L & lin <= length(C)
  if (any(valid)) {
    C[lin[valid]] <- 1
  }

  return(C)
}


# generates an initial matrix with colonized cells according to a layer and points
initial_colonized <- function(data, long, lat, base_layer, proportion = 1, rule = "random",
                              set_seed = 1) {
  l_meta <- layer_metadata(base_layer)

  C <- set_pop(data, long, lat, l_meta$NW_vertex, l_meta$layer_dim, l_meta$cell_size,
               proportion, rule, set_seed)

  return(C)
}



# sample cells from suitable areas
suitable_cells <- function(suit_layer, data = NULL, long = "x", lat = "y") {

  sp <- ifelse(is.null(data), "Species", as.character(data[1, 1]))

  if (is.null(data)) {
    noz <- terra::as.data.frame(suit_layer, xy = TRUE)
    noz <- noz[noz[, 3] > 0, ]

    return(data.frame(species = sp, longitude = noz[, 1], latitude = noz[, 2],
                      suitability = noz[, 3]))
  } else {
    suit_bar <- terra::extract(suit_layer, data[, c(long, lat)])[, 2]
    tokeep <- suit_bar > 0 & !is.na(suit_bar)

    if (sum(tokeep) == 0) {
      stop("None of the records lay in suitable areas to start the simulation.",
           "\n  Try decreasing 'suitability_threshold' or, if projections are done,",
           "\n  try distinct combinations of the periods for the simulation.")
    }

    return(data.frame(species = sp, longitude = data[tokeep, 2],
                      latitude = data[tokeep, 3],
                      suitability = suit_bar[tokeep]))
  }
}


# id cells that have been colonized in a matrix
which_colonized <- function(colonized_matrix, suitability_matrix = NULL,
                            proportion = 1, rule = "random", set_seed = 1) {
  if (!is.null(suitability_matrix)) {
    colonized_cells <- which(colonized_matrix >= 1 & suitability_matrix > 0,
                             arr.ind = TRUE)
  } else {
    colonized_cells <- which(colonized_matrix >= 1, arr.ind = TRUE)
  }

  if (proportion < 1) {
    n <- nrow(colonized_cells)
    ns <- ceiling(n * proportion)

    set.seed(set_seed)
    if (rule == "random") {
      colonized_cells <- colonized_cells[sample(n, ns), ]
    } else {
      suit <- suitability_matrix[colonized_cells]
      colonized_cells <- colonized_cells[sample(n, ns, prob = suit), ]
    }

    if (ns == 1) {
      colonized_cells <- matrix(colonized_cells, ncol = 2)
      colnames(colonized_cells) <- c("row", "col")
    }
  }

  return(colonized_cells)
}


# derive angle and distance of a movement based on a dispersal kernel
angle_distance <- function(n = 1, dispersal_kernel = "normal", kernel_spread = 1,
                           set_seed = 1) {
  set.seed(set_seed)
  theta <- runif(n = n, min = 0, max = 2) * pi

  if (dispersal_kernel == "normal") {
    rad <- rnorm(n = n, mean = 0, sd = kernel_spread)
  }
  if (dispersal_kernel == "log_normal") {
    rad <- rlnorm(n = n, meanlog = 0, sdlog = kernel_spread)
  }

  return(list(angle = theta, distance = rad))
}


# get accessed cells based on angle_distance and previously
which_accessed <- function(colonized_cells, angle_distance, layer_dim) {
  d_lat <- round(angle_distance$distance * sin(angle_distance$angle))
  d_lon <- round(angle_distance$distance * cos(angle_distance$angle))

  r_new <- colonized_cells[, 1] + d_lat
  c_new <- colonized_cells[, 2] + d_lon

  cond <- r_new >= 1 & r_new <= layer_dim[1] & c_new >= 1 & c_new <= layer_dim[2]

  r_final <- colonized_cells[, 1]
  c_final <- colonized_cells[, 2]

  r_final[cond] <- r_new[cond]
  c_final[cond] <- c_new[cond]

  return(cbind(r_final, c_final))
}


# # runs simulation steps to update accessed areas
# dispersal_steps <- function(colonized_matrix, suitability_matrix,
#                             disperser_rules, proportion_to_disperse,
#                             sampling_rule, dispersal_kernel,
#                             kernel_spread, set_seed) {
#
#   # initial tests
#   if (missing(colonized_matrix)) {
#     stop("Argument 'colonized_matrix' must be defined")
#   }
#   if (missing(suitability_matrix)) {
#     stop("Argument 'suitability_matrix' must be defined")
#   }
#   if (missing(disperser_rules)) {
#     stop("Argument 'disperser_rules' must be defined")
#   }
#   if (!dispersal_kernel %in% c("normal", "log_normal")) {
#     stop("Argument 'dispersal_kernel' not valid")
#   }
#   if (!sampling_rule %in% c("random", "suitability")) {
#     stop("Argument 'sampling_rule' is not valid")
#   }
#
#   # # 1. Calcula o número de dispersores para cada célula colonizada
#   # # Nota: Você pode precisar aplicar o nd_sval_opt apenas nos pixels onde C > 0
#   # # Para simplicidade, vamos criar a matriz de dispersores:
#   # Nd_vec <- nd_sval(suitability_matrix, disperser_rules)
#   # Nd_matrix <- matrix(Nd_vec, nrow=nrow(suitability_matrix), ncol=ncol(suitability_matrix))
#
#   # Preparação (o C++ precisa de inputs limpos)
#   is_fat <- (dispersal_kernel == "log_normal")
#
#   # Chama a função C++ que compilamos
#   res <- cpp_dispersal_steps(colonized_matrix, suitability_matrix,
#                              1, kernel_spread, is_fat)
#
#   return(res)
# }

# runs simulation steps to update accessed areas
dispersal_steps <- function(colonized_matrix, suitability_matrix,
                            disperser_rules, proportion_to_disperse,
                            sampling_rule, dispersal_kernel,
                            kernel_spread, set_seed) {

  # initial tests
  if (missing(colonized_matrix)) {
    stop("Argument 'colonized_matrix' must be defined")
  }
  if (missing(suitability_matrix)) {
    stop("Argument 'suitability_matrix' must be defined")
  }
  if (missing(disperser_rules)) {
    stop("Argument 'disperser_rules' must be defined")
  }
  if (!dispersal_kernel %in% c("normal", "log_normal")) {
    stop("Argument 'dispersal_kernel' not valid")
  }
  if (!sampling_rule %in% c("random", "suitability")) {
    stop("Argument 'sampling_rule' is not valid")
  }

  layer_dim <- dim(colonized_matrix)
  A_now <- matrix(0, nrow = layer_dim[1], ncol = layer_dim[2])

  # 1. celulas colonizadas e adequadas (amostragem e seed como no R original)
  colonized_cells <- which_colonized(colonized_matrix, suitability_matrix,
                                     proportion_to_disperse, sampling_rule,
                                     set_seed)
  if (nrow(colonized_cells) == 0) {
    return(A_now)
  }

  # 2. numero de dispersores por celula (1..max_dispersers)
  Sv <- suitability_matrix[colonized_cells]
  n_dispersers <- nd_sval(Sv, disperser_rules)
  Ndmax <- max(n_dispersers)
  total <- sum(n_dispersers)
  if (total == 0) {
    return(A_now)
  }

  # 3. sorteia angulo/distancia no R, rodada a rodada, com o MESMO esquema
  #    de seeds da implementacao original (set_seed + l - 1)
  idx   <- integer(total)
  theta <- numeric(total)
  rad   <- numeric(total)
  pos <- 0L

  for (l in seq_len(Ndmax)) {
    rl <- which(n_dispersers >= l)
    nv <- length(rl)
    if (nv == 0) next

    ad <- angle_distance(n = nv, dispersal_kernel = dispersal_kernel,
                         kernel_spread = kernel_spread,
                         set_seed = set_seed + l - 1)

    w <- pos + seq_len(nv)
    idx[w]   <- rl
    theta[w] <- ad$angle
    rad[w]   <- ad$distance
    pos <- pos + nv
  }

  # C++
  cpp_accumulate_dispersal(A_now, colonized_cells, idx, theta, rad, TRUE)
}

# update A after dispersal step
update_accessed <- function(accessed_matrix, accessed_matrix_now) {
  return(accessed_matrix + accessed_matrix_now)
}


# update C after dispersal step
update_colonized <- function(colonized_matrix, accessed_matrix,
                             suitability_matrix) {
  whichC <- which(accessed_matrix > 0 & suitability_matrix > 0)
  colonized_matrix[whichC] <- colonized_matrix[whichC] + accessed_matrix[whichC]

  return(colonized_matrix)
}
