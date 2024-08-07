#' Simulation of species dispersal processes
#'
#' @description multistep simulation of species dispersal to reconstruct areas
#' that have been or could be accessed and/or colonized based on environmental
#' suitability and user-defined dispersal parameters.
#'
#' @param data data.frame containing geographic coordinates of occurrences of
#' the species of interest. Mandatory columns are: "species", "longitude",
#' "latitude", in that order. Optionally, a fourth column "suitability" should
#' be used when \code{sampling_rule} = "suitability".
#' @param suit_layers (character) vector of names of suitability layers to be
#' used as distinct scenarios. If more than one, the layer names must be ordered
#' starting from first to last scenario. All layers must have the same extent,
#' resolution, number of cells, and projection. Layer names should include
#' parent directory if needed.
#' @param starting_proportion (numeric) proportion of \code{data} to be used as
#' starting points for the simulation. Default = 0.5. All data is used if a
#' value of 1 is defined.
#' @param proportion_to_disperse (numeric) proportion of colonized cells from
#' which dispersers will start a new dispersal event; default = 1.
#' @param sampling_rule (character) rule to be used to sample a
#' \code{starting_proportion} of \code{data} and a \code{proportion_to_disperse}
#' from colonized cells to run dispersal simulation steps. Options are: "random"
#' and "suitability". Using the option "suitability" prioritizes records in
#' with higher suitability values. Default = "random".
#' @param dispersal_kernel (character) dispersal kernel (dispersal function)
#' used to simulate the movement of the species. Options are: "normal",
#' "log_normal". Default = "normal".
#' @param kernel_spread (numeric) standard deviation for the
#' \code{dispersal_kernel}. Default = 1.
#' @param max_dispersers (numeric) maximum number of dispersers that depart from
#' each colonized pixel. Depending on suitability this number will automatically
#' decrease in areas with low suitability values. Default = 4.
#' @param dispersal_events (numeric) number of dispersal events to happen per
#' scenario. A vector of multiple values could be used to define different
#' dispersal events for distinct scenarios. See details; default = 25.
#' @param replicates (numeric) number of times to repeat the simulation
#' per scenario or dispersal event, depending on \code{results_by}. Default = 10.
#' @param threshold (numeric) percentage to be considered when excluding
#' accessed or colonized cells with lower values. Default = 5.
#' @param results_by (character) how to run replicates and return results.
#' Options are: "scenario" or "event". If "scenario" the replicates are produced
#' for each of the scenarios represented by layers in \code{suit_layers}. If
#' "event", replicates are produced per dispersal event. Default = "scenario"
#' @param return (character) the results to be returned or written. Options are:
#' "all", "accessed", "colonized". Default = "all"
#' @param set_seed (numeric) a seed to be used when sampling \code{data}
#' according to \code{starting_proportion}. Default = 1.
#' @param write_to_directory (logical) whether to write results in
#' \code{output_directory}. Default = FALSE.
#' @param write_all (logical) valid if \code{results_by} = "scenario" and
#' \code{write_to_directory} = TRUE. Whether or not to write results for all
#' scenarios. The default, FALSE, writes only the final results.
#' @param raster_format (character) format to use for raster layers to be
#' written. Options are: "GTiff" and "ascii". Default = "GTiff.
#' @param output_directory (character) name of the output directory where
#' results should be written. If this directory does not exist, it will be
#' created.
#' @param parallel (logical) whether to run replicates in parallel.
#' Default = FALSE.
#' @param cores (numeric) number of cores to run replicates in parallel. Only
#' works if parallel = TRUE. Default = 4
#' @param overwrite (logical) whether or not to overwrite the
#' \code{output_directory} if it already exists. Default = FALSE.
#' @param progress_bar (logical) whether or not to show progress bar when
#' running replicates in parallel. Default = TRUE.
#'
#' @export
#' @importFrom terra rast ext writeRaster as.matrix res nrow ncol
#' @importFrom stats quantile runif rnorm rlnorm
#' @importFrom utils read.csv write.csv write.table
#' @importFrom foreach foreach `%dopar%`
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#'
#' @rdname dispersal_simulationR
#'
#' @usage
#' dispersal_simulationR(data, suit_layers, starting_proportion = 0.5,
#'                       proportion_to_disperse = 1, sampling_rule = "random",
#'                       dispersal_kernel = "normal",
#'                       kernel_spread = 1, max_dispersers = 4,
#'                       dispersal_events = 25, replicates = 10,
#'                       threshold = 5, results_by = "scenario",
#'                       return = "all", set_seed = 1,
#'                       write_to_directory = FALSE, write_all = FALSE,
#'                       raster_format = "GTiff", output_directory,
#'                       parallel, cores, overwrite, progress_bar)
#'
#' @return
#' If \code{return} = "all', all elements described below will be returned as a
#' list. If "accessed" or "colonized" are chosen instead, only the elements
#' corresponding to either "accessed" or "colonized" areas will be returned.
#'
#' Notice that only final results are returned as part of the list. If views of
#' all scenarios are needed use options to write results in a directory,
#' \code{write_to_directory}, \code{write_all}, \code{raster_format}, and
#' \code{output_directory}.
#'
#' The list returned contains:
#' - Summary: a list with a summary of scenarios and parameters used
#' - if \code{results_by} = "scenario":
#'     - A: a binary SpatRaster showing accessed = 1 and non-accessed areas = 0
#'     - A_mean: a SpatRaster representing mean values of accessibility
#'     frequency among replicates
#'     - A_var: a SpatRaster representing variance among values of accessibility
#'     frequency of all replicates
#'     - A_scenarios: a SpatRaster with values representing the number of the
#'     scenario when areas where accessed
#'     - C: a binary SpatRaster showing colonized = 1 and non-colonized areas = 0
#'     - C_mean: a SpatRaster representing mean values of frequency of
#'     colonization among replicates
#'     - C_var: a SpatRaster representing variance among values of frequency of
#'     colonization of all replicates
#'     - C_scenarios: a SpatRaster with values representing the number of the
#'     scenario when areas where colonized
#' - if \code{results_by} = "event":
#'     - A_events: a SpatRaster with values representing the number of the
#'     dispersal event when areas where accessed
#'     - C_events: a SpatRaster with values representing the number of the
#'     dispersal event when areas where colonized
#'
#' The number of dispersal events in results is continuous among scenarios. If
#' 10 dispersal events are defined and multiple scenarios exist in
#' \code{suit_layers}, the first dispersal event in the second scenario will be
#' number 11.
#'
#' If \code{write_to_directory} is set to TRUE, the elements described above
#' and raster layers corresponding to all scenarios (if \code{write_all} = TRUE),
#' or events per scenario (depending on \code{results_by}) will be written in
#' \code{output_directory}.
#'
#' @details
#' Defining a vector of multiple values in \code{dispersal_events} could be
#' useful when distinct scenarios represent different periods of time, or if
#' a reduced number of events need to be simulated in the last scenario.
#' If a vector of values is defined in \code{dispersal_events}, the length of
#' this vector must match the length of \code{suit_layers}, otherwise,
#' the first element in \code{dispersal_events} will be used and a
#' message will be printed.
#'
#' @examples
#' # data
#' data("records", package = "grinnell")
#' suitability <- system.file("extdata/suitability.tif", package = "grinnell")
#'
#' # simulation
#' d_s <- dispersal_simulationR(data = records, suit_layers = suitability,
#'                              replicates = 3, dispersal_events = 5,
#'                              return = "colonized")


dispersal_simulationR <- function(data, suit_layers, starting_proportion = 0.5,
                                  proportion_to_disperse = 1,
                                  sampling_rule = "random",
                                  dispersal_kernel = "normal",
                                  kernel_spread = 1, max_dispersers = 4,
                                  dispersal_events = 25, replicates = 10,
                                  threshold = 5, results_by = "scenario",
                                  return = "all", set_seed = 1,
                                  write_to_directory = FALSE, write_all = FALSE,
                                  raster_format = "GTiff", output_directory,
                                  parallel, cores, overwrite, progress_bar) {

  # initial tests
  if (missing(data)) {
    stop("Argument 'data' must be defined")
  }
  if (missing(suit_layers)) {
    stop("Argument 'suit_layers' must be defined")
  }
  if (write_to_directory & missing(output_directory)) {
    stop("If 'write_to_directory' = TRUE, 'output_directory' must be defined")
  }
  if (!sampling_rule %in% c("random", "suitability")) {
    stop("Argument 'sampling_rule' is not valid")
  } else {
    if (sampling_rule == "suitability" & ncol(data) < 4) {
      stop("A fourth column in 'data' containing suitability values is required")
    }
  }
  if (length(dispersal_events) == 1) {
    dispersal_events <- rep(dispersal_events, length(suit_layers))
  } else {
    if (length(dispersal_events) != length(suit_layers)) {
      message("Length of 'dispersal_events' does not match number of scenarios",
              " using first value: ", dispersal_events[1])
      dispersal_events <- rep(dispersal_events[1], length(suit_layers))
    }
  }


  # initial part of report
  if (write_to_directory == TRUE) {
    if (dir.exists(output_directory) == FALSE) {
      dir.create(output_directory)
    }
    outText <- paste0(output_directory, "/report.txt")
    if (file.exists(outText)) {invisible(file.remove(outText))}
    cat(
      "Simulation parameters\n\n",
      "   Suitability scenarios: ", length(suit_layers), "\n",
      "   Dispersal events: ", dispersal_events, "\n",
      "   Replicates: ", replicates, "\n",
      "   Dispersal kernel: ", dispersal_kernel, "\n",
      "   Kernel spread (SD): ", kernel_spread, "\n",
      "   Maximum number of dispersers: ", max_dispersers, "\n",
      file = outText,
      append = TRUE
    )
  }

  # running
  ## start time
  start <- Sys.time()

  ## running all analysis
  message("Running simulation:")

  if (results_by == "scenario") {
    res <- scenario_wise_simulation(data, suit_layers, starting_proportion,
                                    proportion_to_disperse, sampling_rule,
                                    dispersal_kernel, kernel_spread,
                                    max_dispersers, dispersal_events, replicates,
                                    threshold, return, set_seed,
                                    write_to_directory, write_all,
                                    raster_format, output_directory,
                                    cores = cores, parallel,
                                    overwrite = overwrite, progress_bar)
  } else {
    res <- event_wise_simulation(data, suit_layers, starting_proportion,
                                 proportion_to_disperse, sampling_rule,
                                 dispersal_kernel, kernel_spread, max_dispersers,
                                 dispersal_events, replicates,
                                 threshold, return, set_seed,
                                 write_to_directory, raster_format,
                                 output_directory, cores = cores,
                                 parallel, overwrite = overwrite, progress_bar)
  }

  # end time
  end <- Sys.time()

  # last parts of report and preparing layers if needed
  if (write_to_directory == TRUE) {
    timetot <- end - start
    cat("\nSimulation time\n\n",
        "   Start date | time: ", format(start, usetz = TRUE), "\n",
        "   Running time: ", timetot, attr(timetot, "units"),
        file = outText, append = TRUE)
  }

  # preparing results
  res <- list(Summary = res$Summary, A = res$A, A_mean = res$A_mean,
              A_var = res$A_var, A_scenarios = res$A_scenarios,
              A_events = res$A_events, C = res$C, C_mean = res$C_mean,
              C_var = res$C_var, C_scenarios = res$C_scenarios,
              C_events = res$C_events)

  # returning results
  return(res)
}



#' @rdname dispersal_simulationR
#' @export
#' @usage
#' scenario_wise_simulation(data, suit_layers, starting_proportion = 0.5,
#'                          proportion_to_disperse = 1, sampling_rule = "random",
#'                          dispersal_kernel = "normal",
#'                          kernel_spread = 1, max_dispersers = 4,
#'                          dispersal_events = 25, replicates = 10,
#'                          threshold = 5, return = "all",
#'                          set_seed = 1, write_to_directory = FALSE,
#'                          write_all = FALSE, raster_format = "GTiff",
#'                          output_directory)

scenario_wise_simulation <- function(data, suit_layers, starting_proportion = 0.5,
                                     proportion_to_disperse = 1,
                                     sampling_rule = "random",
                                     dispersal_kernel = "normal",
                                     kernel_spread = 1, max_dispersers = 4,
                                     dispersal_events = 25, replicates = 10,
                                     threshold = 5, return = "all",
                                     set_seed = 1, write_to_directory = FALSE,
                                     write_all = FALSE, raster_format = "GTiff",
                                     output_directory, parallel = TRUE,
                                     cores, overwrite, progress_bar) {

  # initial values
  cur_layer <- terra::rast(suit_layers[length(suit_layers)])
  l_meta <- layer_metadata(cur_layer)
  layer_dim <- l_meta$layer_dim

  # other parameters
  d_rules <- disperser_rules(cur_layer, max_dispersers)

  # layers to store when things happen
  ns <- length(suit_layers) + 1
  cur_layer[] <- ns

  if(return == "all") {
    a_when <- cur_layer; c_when <- cur_layer
  } else {
    if(return == "accessed") {
      a_when <- cur_layer
    } else {
      c_when <- cur_layer
    }
  }

  # lists of accessed and colonized
  list_acc <- list()
  list_col <- list()

  # simulation
  for (i in 1:length(suit_layers)) {
    message("  Scenario ", i, " of ", length(suit_layers),
            appendLF = parallel) #Progress bar depends on parallel

    s <- terra::rast(suit_layers[i])
    S <- base_matrix(s)

    #If not in parallel (or 1 core or only 1 replicate)
    if(isFALSE(parallel) | cores == 1 | replicates == 1) {

    ## loop for all replicates
    message(" - Replicate:", appendLF = FALSE)

    for (j in 1:replicates) {
      ## preparing C matrix
      set_seed <- set_seed + j - 1
      if (i == 1) {
        C <- set_pop(data, long = "longitude", lat = "latitude",
                     l_meta$NW_vertex, layer_dim, l_meta$cell_size,
                     starting_proportion, sampling_rule, set_seed)
        A <- C
      } else {
        A <- matrix(list_acc[[j]], nrow = layer_dim[1], ncol = layer_dim[2])
        C <- matrix(list_col[[j]], nrow = layer_dim[1], ncol = layer_dim[2])
        C[S == 0] <- 0
      }

      ### simulation steps
      for (k in 1:dispersal_events[i]) {
        #### running steps according to dispersers
        set_seed1 <- set_seed + ((k - 1) * 10)
        A_now <- dispersal_steps(colonized_matrix = C, suitability_matrix = S,
                                 disperser_rules = d_rules, proportion_to_disperse,
                                 sampling_rule, dispersal_kernel, kernel_spread,
                                 set_seed1)

        #### updating A and C
        A <- update_accessed(A, A_now)
        C <- update_colonized(C, A_now, S)
      }

      ### keeping replicates
      list_acc[[j]] <- c(A)
      list_col[[j]] <- c(C)
      message(" ", j, appendLF = FALSE)

      }} else { #In parallel

        #Make ans register cluster of cores
        cl <- parallel::makeCluster(cores)
        doSNOW::registerDoSNOW(cl)
        `%dopar%` <- foreach::`%dopar%` # so %dopar% doesn't need to be attached

        #Show progress bar?
        if (progress_bar) {
          pb <- txtProgressBar(min = 0, max = replicates, style = 3)
          progress <- function(n) setTxtProgressBar(pb, n)
          opts <- list(progress = progress)} else {
            opts <- NULL
          }

        result_replicates <- foreach::foreach(j = 1:replicates,
            .options.snow = opts,
            .export = c("set_pop", "set_loc", "dispersal_steps",
                        "which_colonized", "nd_sval", "angle_distance",
                        "which_accessed", "update_accessed",
                        "update_colonized", "replicate_stats",
                        "binarize_matrix")) %dopar% {
      set_seed <- set_seed + j - 1
      if (i == 1) {
        C <- set_pop(data, l_meta$NW_vertex, layer_dim,
                     l_meta$cell_size, starting_proportion, sampling_rule,
                     set_seed)
        A <- C
      }      else {
        A <- matrix(list_acc[[j]], nrow = layer_dim[1],
                    ncol = layer_dim[2])
        C <- matrix(list_col[[j]], nrow = layer_dim[1],
                    ncol = layer_dim[2])
        C[S == 0] <- 0
      }
      for (k in 1:dispersal_events[i]) {
        set_seed1 <- set_seed + ((k - 1) * 10)
        A_now <- dispersal_steps(colonized_matrix = C,
                                 suitability_matrix = S, disperser_rules = d_rules,
                                 proportion_to_disperse, sampling_rule, dispersal_kernel,
                                 kernel_spread, set_seed1)
        A <- update_accessed(A, A_now)
        C <- update_colonized(C, A_now, S)
      }
      return(list("A" = as.numeric(A), "C" = as.numeric(C)))
    }

    list_acc <- lapply(result_replicates, function(x) x$A)
    list_col <- lapply(result_replicates, function(x) x$C)
    parallel::stopCluster(cl)
    }

    ### statistics an updates and correcting with suitability
    s <- (s > 0) * 1
    if(return == "all") {
      Amvb <- replicate_stats(list_acc, S, s, threshold)
      Cmvb <- replicate_stats(list_col, S, s, threshold)

      a_when <- a_when - Amvb[[3]]
      c_when <- (c_when - Cmvb[[3]]) * s
    } else {
      if(return == "accessed") {
        Amvb <- replicate_stats(list_acc, S, s, threshold)
        a_when <- a_when - Amvb[[3]]
      } else {
        Cmvb <- replicate_stats(list_col, S, s, threshold)
        c_when <- (c_when - Cmvb[[3]]) * s
      }
    }

    ## writing results
    if (write_to_directory == TRUE) {
      if (write_all == TRUE) {
        if(return == "all") {
          namesA <- paste0("A_", c("mean", "var", "bin"), "_scenario_", i)
          write_stats(Amvb, namesA, raster_format, output_directory, overwrite)
          namesC <- paste0("C_", c("mean", "var", "bin"), "_scenario_", i)
          write_stats(Cmvb, namesC, raster_format, output_directory, overwrite)
        } else {
          if(return == "accessed") {
            namesA <- paste0("A_", c("mean", "var", "bin"), "_scenario_", i)
            write_stats(Amvb, namesA, raster_format, output_directory,
                        overwrite)
          } else {
            namesC <- paste0("C_", c("mean", "var", "bin"), "_scenario_", i)
            write_stats(Cmvb, namesC, raster_format, output_directory,
                        overwrite)
          }
        }
      } else {
        if(i == length(suit_layers)) {
          if(return == "all") {
            namesA <- paste0("A_", c("mean", "var", "bin"))
            write_stats(Amvb, namesA, raster_format, output_directory,
                        overwrite)
            namesC <- paste0("C_", c("mean", "var", "bin"))
            write_stats(Cmvb, namesC, raster_format, output_directory,
                        overwrite)
          } else {
            if(return == "accessed") {
              namesA <- paste0("A_", c("mean", "var", "bin"))
              write_stats(Amvb, namesA, raster_format, output_directory,
                          overwrite)
            } else {
              namesC <- paste0("C_", c("mean", "var", "bin"))
              write_stats(Cmvb, namesC, raster_format, output_directory,
                          overwrite)
            }
          }
        }
      }
    }

    message(" #")
  }

  # returning results
  form1 <- rformat_type(raster_format)

  summ <- list(Scenarios = length(suit_layers),
               starting_proportion = starting_proportion,
               proportion_to_disperse = proportion_to_disperse,
               sampling_rule = sampling_rule,
               Dispersal_events = dispersal_events, Replicates = replicates,
               Dispersal_kernel = dispersal_kernel,
               Kernel_spread_SD = kernel_spread, Max_dispersers = max_dispersers)

  if(return == "all") {
    a_when[a_when[] == ns] <- 0
    c_when[c_when[] == ns] <- 0
    c_when[c_when[] < 0] <- 0

    if (write_to_directory == TRUE) {
      aname <- paste0(output_directory, "/A_classified", form1)
      cname <- paste0(output_directory, "/C_classified", form1)

      terra::writeRaster(a_when, filename = aname, overwrite = overwrite)
      terra::writeRaster(c_when, filename = cname, overwrite = overwrite)
    }

    res <- list(Summary = summ, A = Amvb[[3]], A_mean = Amvb[[1]],
                A_var = Amvb[[2]], A_scenarios = a_when, C = Cmvb[[3]],
                C_mean = Cmvb[[1]], C_var = Cmvb[[2]], C_scenarios = c_when)
  } else {
    if(return == "accessed") {
      a_when[a_when[] == ns] <- 0

      if (write_to_directory == TRUE) {
        aname <- paste0(output_directory, "/A_scenarios", form1)
        terra::writeRaster(a_when, filename = aname, overwrite = overwrite)
      }

      res <- list(Summary = summ, A = Amvb[[3]], A_mean = Amvb[[1]],
                  A_var = Amvb[[2]], A_scenarios = a_when, C = NULL,
                  C_mean = NULL, C_var = NULL, C_scenarios = NULL)
    } else {
      c_when[c_when[] == ns] <- 0
      c_when[c_when[] < 0] <- 0

      if (write_to_directory == TRUE) {
        cname <- paste0(output_directory, "/C_scenarios", form1)
        terra::writeRaster(c_when, filename = cname, overwrite = overwrite)
      }

      res <- list(Summary = summ, A = NULL, A_mean = NULL, A_var = NULL,
                  A_scenarios = NULL, C = Cmvb[[3]], C_mean = Cmvb[[1]],
                  C_var = Cmvb[[2]], C_scenarios = c_when)
    }
  }

  return(res)
}



#' @rdname dispersal_simulationR
#' @export
#' @usage
#' event_wise_simulation(data, suit_layers, starting_proportion = 0.5,
#'                       proportion_to_disperse = 1, sampling_rule = "random",
#'                       dispersal_kernel = "normal",
#'                       kernel_spread = 1, max_dispersers = 4,
#'                       dispersal_events = 25, replicates = 10,
#'                       threshold = 5, return = "all",
#'                       set_seed = 1, write_to_directory = FALSE,
#'                       raster_format = "GTiff", output_directory)

event_wise_simulation <- function(data, suit_layers, starting_proportion = 0.5,
                                  proportion_to_disperse = 1,
                                  sampling_rule = "random",
                                  dispersal_kernel = "normal",
                                  kernel_spread = 1, max_dispersers = 4,
                                  dispersal_events = 25, replicates = 10,
                                  threshold = 5, return = "all",
                                  set_seed = 1, write_to_directory = FALSE,
                                  raster_format = "GTiff", output_directory,
                                  overwrite = TRUE) {

  # initial values
  cur_layer <- terra::rast(suit_layers[length(suit_layers)])
  l_meta <- layer_metadata(cur_layer)
  layer_dim <- l_meta$layer_dim

  # other parameters
  d_rules <- disperser_rules(cur_layer, max_dispersers)

  # layers to store when things happen
  ne <- sum(dispersal_events) + 1
  cur_layer[] <- ne

  if(return == "all") {
    a_when <- cur_layer; c_when <- cur_layer
  } else {
    if(return == "accessed") {
      a_when <- cur_layer
    } else {
      c_when <- cur_layer
    }
  }

  # lists of accessed and colonized
  list_acc <- list()
  list_col <- list()

  # simulation
  for (i in 1:length(suit_layers)) {
    message("  Scenario ", i, " of ", length(suit_layers), appendLF = FALSE)

    s <- terra::rast(suit_layers[i])
    S <- base_matrix(s)

    ## loop for all steps
    message(" - D. event:", appendLF = FALSE)

    for (j in 1:dispersal_events[i]) {

      for (k in 1:replicates) {
        ### preparing C matrix
        set_seed <- set_seed + k - 1
        if (i == 1 & j == 1) {
          C <- set_pop(data, l_meta$NW_vertex, layer_dim, l_meta$cell_size,
                       starting_proportion, sampling_rule, set_seed)
          A <- C
        } else {
          A <- matrix(list_acc[[k]], nrow = layer_dim[1], ncol = layer_dim[2])
          C <- matrix(list_col[[k]], nrow = layer_dim[1], ncol = layer_dim[2])
          C[S == 0] <- 0
        }

        ### updating matrices
        #### running steps according to dispersers
        A_now <- dispersal_steps(colonized_matrix = C, suitability_matrix = S,
                                 disperser_rules = d_rules, proportion_to_disperse,
                                 sampling_rule, dispersal_kernel, kernel_spread,
                                 set_seed)

        #### updating A and C
        A <- update_accessed(A, A_now)
        C <- update_colonized(C, A_now, S)

        ### keeping replicates
        list_acc[[k]] <- c(A)
        list_col[[k]] <- c(C)

      }

      ### statistics an updates
      if(return == "all") {
        Amvb <- replicate_stats(list_acc, S, s, threshold)
        Cmvb <- replicate_stats(list_col, S, s, threshold)

        a_when <- a_when - Amvb[[3]]
        c_when <- c_when - Cmvb[[3]]
      } else {
        if(return == "accessed") {
          Amvb <- replicate_stats(list_acc, S, s, threshold)
          a_when <- a_when - Amvb[[3]]
        } else {
          Cmvb <- replicate_stats(list_col, S, s, threshold)
          c_when <- c_when - Cmvb[[3]]
        }
      }

      ### writing results
      if (write_to_directory == TRUE) {
        if(return == "all") {
          namesA <- paste0("A_", c("mean", "var", "bin"), "_scenario_", i,
                           "_event_", j)
          write_stats(Amvb, namesA, raster_format, output_directory, overwrite)
          namesC <- paste0("C_", c("mean", "var", "bin"), "_scenario_", i,
                           "_event_", j)
          write_stats(Cmvb, namesC, raster_format, output_directory, overwrite)
        } else {
          if(return == "accessed") {
            namesA <- paste0("A_", c("mean", "var", "bin"), "_scenario_", i,
                             "_event_", j)
            write_stats(Amvb, namesA, raster_format, output_directory,
                        overwrite)
          } else {
            namesC <- paste0("C_", c("mean", "var", "bin"), "_scenario_", i,
                             "_event_", j)
            write_stats(Cmvb, namesC, raster_format, output_directory,
                        overwrite)
          }
        }
      }

      message(" ", j, appendLF = FALSE)
    }

    # correcting with suitability
    s <- (s > 0) * 1
    if(return == "all") {
      c_when <- c_when * s
    } else {
      if(return == "colonized") {
        c_when <- c_when * s
      }
    }

    message(" #")
  }

  # returning results
  form1 <- rformat_type(raster_format)

  summ <- list(Scenarios = length(suit_layers),
               starting_proportion = starting_proportion,
               proportion_to_disperse = proportion_to_disperse,
               sampling_rule = sampling_rule,
               Dispersal_events = dispersal_events, Replicates = replicates,
               Dispersal_kernel = dispersal_kernel,
               Kernel_spread_SD = kernel_spread, Max_dispersers = max_dispersers)

  if(return == "all") {
    a_when[a_when[] == ne] <- 0
    c_when[c_when[] == ne] <- 0
    c_when[c_when[] < 0] <- 0

    if (write_to_directory == TRUE) {
      aname <- paste0(output_directory, "/A_classified", form1)
      cname <- paste0(output_directory, "/C_classified", form1)

      terra::writeRaster(a_when, filename = aname, overwrite = overwrite)
      terra::writeRaster(c_when, filename = cname, overwrite = overwrite)
    }

    res <- list(Summary = summ, A_events = a_when, C_events = c_when)
  } else {
    if(return == "accessed") {
      a_when[a_when[] == ne] <- 0

      if (write_to_directory == TRUE) {
        aname <- paste0(output_directory, "/A_events", form1)
        terra::writeRaster(a_when, filename = aname, overwrite = overwrite)
      }

      res <- list(Summary = summ, A_events = a_when, C_events = NULL)
    } else {
      c_when[c_when[] == ne] <- 0
      c_when[c_when[] < 0] <- 0

      if (write_to_directory == TRUE) {
        cname <- paste0(output_directory, "/C_events", form1)
        terra::writeRaster(c_when, filename = cname, overwrite = overwrite)
      }

      res <- list(Summary = summ, A_events = NULL, C_events = c_when)
    }
  }

  return(res)
}
