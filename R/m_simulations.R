#' Simulation of species accessible areas (M)
#' @description
#' m_simulation generates an area that has potentially been accessible to a
#' species based on simulations of dispersal events determined by environmental
#' suitability and user-defined parameters. This version differs from the
#' original M_simulationR as it uses defined past scenarios instead of
#' interpolations, does not involve reading or writing files, includes
#' additional steps to prepare the final M, and is able to test several
#' combinations of dispersal events and kernel spread at once.
#' @param data (data.frame) occurrence records for the species of interest to be
#' used to run the simulation. Columns must include species name, longitude,
#' and latitude
#' @param long (character) column name in data with longitude
#' @param lat (character) column name in data with latitude
#' @param current_variables (SpatRaster) environmental variables representing
#' "current" conditions.
#' @param variables (character) variable names to subset current variables. If
#' NULL (default), use all variables.
#' @param scale (logical) whether or not to scale variables while performing
#' principal component analyses.
#' Default = TRUE.
#' @param center (logical) whether or not to center variables while performing
#' principal component analyses.
#' Default = TRUE.
#' @param projection_dir (character) directory with raster files representing
#' past conditions.
#' @param periods (character) names of the periods representing past conditions,
#' order from the most recent to oldest.
#' @param pattern (character) format of the raster files representing past
#' conditions. Default = ".tif".
#' @param initial_m_buffer (numeric) width (in kilometers) of the buffer to be
#' added to the initial M (See details). If NULL (default), the initial M will
#' have exactly the same extent as the current variables.
#' @param suitability_threshold (numeric) value (in percentage) to be used as
#' threshold for suitability. Below this value environments are considered
#' unsuitable. Default = 5.
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
#' @param comb_grid (data.frame) A grid containing combinations of dispersal events
#' and kernel spreads for running the simulations. The columns must be named
#' "dispersal events" and "kernel_spread". See details for more information.
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
#' @param verbose (logical) whether to display messages.
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
#'     - information on the paremeters used: number of scenarios,
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
#' If multiple combinations of kernel_spread and dispersal_events are used,
#' the function returns a list containing results for each combination.
#' @export
#' @importFrom terra extract
m_simulations <- function(data, long, lat, current_variables, variables,
                         scale = TRUE, center = TRUE,
                         projection_dir, periods, pattern = ".tif",
                         initial_m_buffer = NULL,
                         suitability_threshold = 5,
                         starting_proportion = 0.75,
                         proportion_to_disperse = 1,
                         sampling_rule = "random",
                         dispersal_kernel = "normal",
                         kernel_spread = 2, max_dispersers = 4,
                         dispersal_events = 25,
                         comb_grid = NULL,
                         replicates = 10,
                         threshold = 5,
                         set_seed = 1,
                         skip_extinction = TRUE,
                         results_by_event = FALSE,
                         results_by_scenario = FALSE,
                         remove_m_without_records = TRUE,
                         extra_buffer = 25,
                         progress_bar = TRUE,
                         verbose = TRUE){
  if(verbose){
    message("Predicting suitability models")
  }
  #Get models
  model <- get_models(data, long, lat,
                      current_variables,
                      variables,
                      scale,
                      center,
                      projection_dir,
                      periods,
                      pattern,
                      initial_m_buffer,
                      suitability_threshold)

  #Keep one record by pixel
  new_occ <- data
  new_occ$cell <- terra::extract(model$models[[1]], data[,c(long, lat)],
                          cell = TRUE, ID = FALSE)$cell
  new_occ <- new_occ[!duplicated(new_occ$cell), ]
  new_occ$cell <- NULL

  if(verbose){
    message("Simulating dispersion in past scenarios...")
  }

  if(is.null(comb_grid)){
    comb_grid <- data.frame(dispersal_events = dispersal_events,
                            kernel_spread = kernel_spread)
  }


  #Create list to save simulations
  res <- list()
  for(i in 1:nrow(comb_grid)){
    #print(i)
    if(verbose & nrow(comb_grid) > 1){message("Simulating combination ", i)}

    #Simulation
    ms <- scenario_wise_simulation_b(data, long, lat,
                                     suit_layers = model$models,
                                     starting_proportion,
                                     proportion_to_disperse,
                                     sampling_rule,
                                     dispersal_kernel,
                                     kernel_spread = comb_grid[i, "kernel_spread"],
                                     max_dispersers,
                                     dispersal_events = comb_grid[i, "dispersal_events"],
                                     replicates,
                                     threshold,
                                     set_seed,
                                     skip_extinction,
                                     results_by_event,
                                     results_by_scenario,
                                     remove_m_without_records,
                                     extra_buffer,
                                     progress_bar)
    res[[i]] <- ms
  }
  names(res) <- paste0("Combination_", 1:nrow(comb_grid))
  if(length(res) == 1){
    res <- res[[1]]
  } else{
    all_summary <- do.call("rbind", lapply(res, function(x) x$summary))
    row.names(all_summary) <- NULL
    for(x in 1:length(res)){
      res[[x]]$summary <- NULL}
    res$summary <- all_summary
  }
  #Append suitability models
  res$suitability_models <- model$models

  return(res)
}
