#' Simulation of species dispersal processes - Version B
#'
#' @description Multistep simulation of species dispersal to reconstruct areas
#' that have been or could be accessed and/or colonized based on environmental
#' suitability and user-defined dispersal parameters. This version differs from
#' the original dispersal_simulation as it does not involve reading or writing
#' files and includes additional steps to prepare the final M (see
#' `remove_m_without_records` and `extra_buffer` parameters).
#'
#' @param data (data.frame) occurrence records for the species of interest to be
#' used to run the simulation. Columns must include species name, longitude,
#' and latitude
#' @param long (character) column name in data with longitude
#' @param lat (character) column name in data with latitude
#' @param suit_layers (SpatRaster) suitability layers to be used as distinct
#' scenarios. If more than one, the layer names must be combined in a single
#' spatraster object, ordered starting from current to last past scenario.
#' @param starting_proportion (numeric) proportion of data to be used as
#' starting points for the simulation. Default = 0.5. All data is used if a
#' value of 1 is defined.
#' @param proportion_to_disperse (numeric) proportion of colonized cells from
#' which dispersers will start a new dispersal event; default = 1.
#' @param sampling_rule (character) rule to be used to sample a
#' starting_proportion of data and a proportion_to_disperse from colonized cells
#' to run dispersal simulation steps. Options are: "random" and "suitability".
#' Using the option "suitability" prioritizes records in with higher
#' suitability values. Default = "random".
#' @param dispersal_kernel (character) dispersal kernel (dispersal function)
#' used to simulate the movement of the species. Options are: "normal",
#' "log_normal". Default = "normal".
#' @param kernel_spread (numeric) standard deviation for the dispersal_kernel.
#' Default = 1.
#' @param max_dispersers (numeric) maximum number of dispersers that depart from
#' each colonized pixel. Depending on suitability this number will automatically
#' decrease in areas with low suitability values. Default = 4.
#' @param dispersal_events (numeric) number of dispersal events to happen per
#' scenario. A vector of multiple values could be used to define different
#' dispersal events for distinct scenarios. See details. Default = 25.
#' @param replicates (numeric) number of times to repeat the simulation per
#' scenario. Default = 10
#' @param threshold (numeric) percentage to be considered when excluding
#' accessed or colonized cells with lower values. Default = 5.
#' @param set_seed (numeric) a seed to be used when sampling data according to
#' starting_proportion. Default = 1.
#' @param skip_extinction (logical) whether to skip scenarios where all
#' populations are extinct (i.e., there are no suitable cells able to be colonized).
#' @param results_by_event (logical) wheter to return Spatrasters representing
#' the acessed, colonized and suitable-but-inaccessible cells by event.
#' @param results_by_scenario (logical) wheter to return Spatrasters representing
#' the acessed, colonized and suitable-but-inaccessible cells by scenario.
#' @param remove_m_without_records (logical) whether to remove disjunct polygons
#' that do not have any records inside them.
#' @param extra_buffer (numeric) width (in kilometers) of the buffer to be added
#' to the records outside the M. Only the buffers that connect to the M are kept
#' in the final M. If NULL, no buffers are drawn around outside records.
#' @param progress_bar (logical) whether to show progress bar. Default = TRUE.
#'
#' @return
#' A list containing:
#' - m_final: a SpatVector of the final M, with extra buffer around outside
#' records (if extra_buffer = TRUE) and without disjunct polygons without
#' records (if remove_m_without_records = TRUE)
#' - m_raw: A SpatVector of the raw M, without any extra buffer and including
#' all disjoint polygons.
#' - m_by_event: a combined Spatraster representing the acessed, colonized and
#' suitable-but-inaccessible cells by event. Only returned if
#' results_by_event = TRUE.
#' - m_by_scen: a combined Spatraster representing the acessed, colonized and
#' suitable-but-inaccessible cells by scenario. Only returned if
#' results_by_scenario = TRUE.
#' - summary: a data.frame including:
#'     - information on the paremeters used: nunber of scenarios,
#'     starting_proportion, proportion_to_disperse, sampling_rule,
#'     Dispersal_events, Replicates, Dispersal_kernel, Kernel_spread_SD, and
#'     Max_dispersers.
#'     - n_polygons_with_occ: number of disjoint polygons where the species
#'     occurs within the initial M. For example, if a species occurs in the
#'     mainland and 3 other islands, this value would be 4.
#'     - n_m: Number of disjoint polygons representing the final M.
#'     If n_m > n_polygons_with_occ, it may indicate that the final M consists
#'     of disjoint polygons where there should be joined polygons.
#' - data: The original data records, with an additional column indicating
#'     if the record falls outside the final M.
#'
#' @export
#' @importFrom terra minmax extract nlyr disagg as.polygons trim is.related vect aggregate union expanse crs
#' @importFrom progress progress_bar
#'
scenario_wise_simulation_b <- function(data, long, lat, suit_layers,
                                       starting_proportion = 0.75,
                                       proportion_to_disperse = 1,
                                       sampling_rule = "random",
                                       dispersal_kernel = "normal",
                                       kernel_spread = 2, max_dispersers = 4,
                                       dispersal_events = 25, replicates = 10,
                                       threshold = 5,
                                       set_seed = 1,
                                       skip_extinction = TRUE,
                                       results_by_event = FALSE,
                                       results_by_scenario = FALSE,
                                       remove_m_without_records = TRUE,
                                       extra_buffer = 25,
                                       progress_bar = TRUE) {

  # initial values
  cur_layer <- suit_layers[[1]]
  l_meta <- layer_metadata(cur_layer)
  layer_dim <- l_meta$layer_dim

  #Remove scenarios with extinction
  if(skip_extinction){
    #Get maximum values
    max_val <- sapply(suit_layers, function(x) terra::minmax(x)[2])
    suit_layers <- suit_layers[[max_val > 0]]
  }

  #If sampling rule = "suitability"
  if(sampling_rule == "suitability"){
    data$suitability <- terra::extract(cur_layer, data[,c(long, lat)], ID = FALSE)}

  # other parameters
  d_rules <- disperser_rules(cur_layer, max_dispersers)

  # layers to store when things happen
  ns <- length(suit_layers) + 1
  cur_layer[] <- ns

  # a_when <- cur_layer
  # c_when <- cur_layer

  #### lists of accessed and colonized ####
  list_acc <- list()
  list_col <- list()
  #If return by scenario#
  if(results_by_scenario){
    list_scen <- list()
  } else {list_scen <- NULL}
  #If return by event
  if(results_by_event){
    list_event <- rep(list(NULL), terra::nlyr(suit_layers))
    list_event <- lapply(list_event, function(x) vector("list", dispersal_events))
  } else {list_event <-  NULL}

  if (progress_bar) {
    # Replicates
    total_steps <- replicates * terra::nlyr(suit_layers) * dispersal_events
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent eta::eta [:elapsed]",
      total = total_steps,
      clear = FALSE
    )
  }

  #### simulation ####
  for (i in 1:terra::nlyr(suit_layers)) {
    s <- suit_layers[[i]]
    S <- base_matrix(s)

    #If sampling rule is suitability, get values
    if(sampling_rule == "suitability"){
      suit_values <- terra::extract(cur_layer, data[,c(long, lat)], ID = FALSE,
                                    raw = TRUE)} else {
        suit_values = NULL
      }


    ## loop for all replicates
    for (j in 1:replicates) {
      ## preparing C matrix
      set_seed <- set_seed + i - 1
      if (i == 1) {
        C <- set_pop(data, long, lat, l_meta$NW_vertex, layer_dim, l_meta$cell_size,
                     starting_proportion, sampling_rule, set_seed,
                     suit_values)
        A <- C
      } else {
        A <- matrix(list_acc[[j]], nrow = layer_dim[1], ncol = layer_dim[2])
        C <- matrix(list_col[[j]], nrow = layer_dim[1], ncol = layer_dim[2])
        C[S == 0] <- 0
      }

      ### simulation steps
      for (k in 1:dispersal_events) {
        if (progress_bar) {
          pb$tick()
        }
        #### running steps according to dispersers
        set_seed1 <- set_seed + ((k - 1) * 10)
        #Check nd_sval here: it uses sapply, try to vectorize
        A_now <- try(dispersal_steps(colonized_matrix = C, suitability_matrix = S,
                                                disperser_rules = d_rules, proportion_to_disperse,
                                                sampling_rule, dispersal_kernel, kernel_spread,
                                                set_seed1), silent = TRUE)
        if(is.matrix(A_now)){
          #### updating A and C
          A <- update_accessed(A, A_now)
          C <- update_colonized(C, A_now, S)
        }
        else {
          if(skip_extinction){
            A <- A
            C <- C} else {
              stop("No suitable pixels in scenario ", i, ". Try set skip_extinction = TRUE")
            }}

        if(results_by_event){
          ### statistics an updates and correcting with suitability
          s_event <- (s > 0) * 1
          r_acc <- replicate_stats(list(A), S, s_event, threshold)$bin
          r_acc <- (r_acc >= 1) * 2
          r_col <- replicate_stats(list(C), S, s_event, threshold)$bin
          r_col <- (r_col >= 1) * 3
          r_acc_col <- max(s_event, r_acc, r_col)
          levels(r_acc_col) <- data.frame(id=0:3,
                                          cover=c("Non-accessed",
                                                  "Suitable-but-inaccessible",
                                                  "Accessed", "Colonized"))
          names(r_acc_col) <- paste0("Event ", k, " in ", names(suit_layers)[i])
          list_event[[i]][[k]] <- r_acc_col
        }

      }

      ### keeping replicates
      #### How to save results by event to make animation?
      list_acc[[j]] <- c(A)
      list_col[[j]] <- c(C)
    }

    #If return by scenario#
    if(results_by_scenario){
      ### statistics an updates and correcting with suitability
      s_scen <- (s > 0) * 1
      r_acc <- replicate_stats(list_acc, S,  s_scen, threshold)$bin
      r_acc <- (r_acc >= 1) * 2
      r_col <- replicate_stats(list_col, S,  s_scen, threshold)$bin
      r_col <- (r_col >= 1) * 3
      r_acc_col <- max(s_scen, r_acc, r_col)
      levels(r_acc_col) <- data.frame(id=0:3,
                                      cover=c("Non-accessed",
                                              "Suitable-but-inaccessible",
                                              "Accessed", "Colonized"))
      list_scen[[i]] <- r_acc_col
    }

  }

  ### statistics an updates and correcting with suitability
  s <- (s > 0) * 1

  #Accessed
  Amvb <- replicate_stats(list_acc, S, s, threshold)
  #Colonized
  Cmvb <- replicate_stats(list_col, S, s, threshold)


  #Return final shapefile
    #Get initial M
  initial_m <- terra::disagg(terra::as.polygons(
    terra::trim(!is.na(suit_layers[[1]]) * 1,0)))
  #Crop layer to initial M
  A_bin <- terra::crop(Amvb[[3]], initial_m, mask = TRUE)
  A_bin[A_bin[] == 0] <- NA
  A_bin <- terra::trim(A_bin)
  shpm <- terra::as.polygons(A_bin, dissolve = TRUE)
  #plot(shpm)
  #mapview(shpm)

  if(results_by_scenario){
    list_scen <- terra::rast(list_scen)
    names(list_scen) <- names(suit_layers)
  }

  if(results_by_event){
    list_event <- terra::rast(unlist(list_event))
  }

  ####Prepare m final ####
  #Remove disjunct M without records
  if(remove_m_without_records){
    pts <- terra::vect(data, geom = c(x = long, y = lat),
                       crs = terra::crs(suit_layers))
    #Get only polygons with occurrences
    m_final <- terra::disagg(shpm)
    m_final <- m_final[terra::is.related(m_final, pts, "intersects")]
  } else {m_final <- shpm}



  #Final buffer around records
  pts_outside <- pts[!terra::is.related(pts, m_final, "intersects")]

  #Add final buffer
  if(!is.null(extra_buffer) & length(pts_outside) > 0) {
    b_occ <- terra::buffer(pts_outside , width = extra_buffer * 1000)
    #Select only buffers that insersects with M
    b_occ <- b_occ[terra::is.related(b_occ, m_final, "intersects")]
    #Merge
    if(length(b_occ) > 0)
    m_final <- terra::aggregate(terra::union(m_final, b_occ))
    #m_v <- buffer(m_v, width = final_buffer * 1000)
  }

  #Get final report
  #Flag records outside the M
  data$inside_m <- terra::is.related(pts, m_final, "intersects")
  #Get final report
  #Number of polygons with occurrence in the initial m
  #Get initial_m
  n_pol_occ <- sum(terra::is.related(initial_m, pts, "intersects"))
  #Number of polygons in m final
  n_pol <- length(terra::disagg(m_final))
  #Number of records outside
  n_out <- sum(!terra::is.related(pts, m_final, "intersects"))

  #Area
  m_area <- sum(terra::expanse(m_final,
                               unit = "km"))

  # Summary of results
  summ <- data.frame(Scenarios = terra::nlyr(suit_layers),
                     starting_proportion = starting_proportion,
                     proportion_to_disperse = proportion_to_disperse,
                     sampling_rule = sampling_rule,
                     Dispersal_events = dispersal_events,
                     Replicates = replicates,
                     Dispersal_kernel = dispersal_kernel,
                     Kernel_spread_SD = kernel_spread,
                     Max_dispersers = max_dispersers,
                     n_polygons_with_occ = n_pol_occ,
                     n_m = n_pol,
                     surplus_m = n_pol - n_pol_occ,
                     occ_outside = n_out,
                     area_km = m_area)

  return(list(m_final = m_final,
              m_raw = shpm,
              m_by_event = list_event,
              m_by_scen = list_scen,
              summary = summ,
              occ = data))
}
