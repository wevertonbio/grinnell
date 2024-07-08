#' Interpolate from the Last Glacial Maximum (LGM) to earlier time periods as oscillayers
#' @description
#' oscillayer generates bioclimatic variables from the Last Glacial Maximum (LGM)
#' to 5.4 million years ago with an interval of 10,000 years, based on the approach
#' developed by Gamisch (2019) on Oscillayers. See details and references for more
#' information.
#'
#' @param current_variables (SpatRaster) environmental variables representing
#' "current" conditions.
#' @param lgm (SpatRaster) environmental variables representing the conditions
#' of the Last Glacial Maximum (LGM). It must have the same variables, variable names,
#' and resolution as the current_variables.
#' @param final_period (numeric) The final period (in thousands of years before
#' present) for which interpolation is performed, with intervals of 10,000 years.
#' Valid values range from 30 to 5380 and must be multiples of 10. Default = 120.
#' @param precipitation_var (character) names of the variables representing precipitation
#' (e.g., c("bio_12", "bio_13", etc.)). This is necessary for fixing negative
#' values in the interpolated precipitation variables.
#' @param write_files (logical) whether to write TIF files to disk. Default = FALSE.
#' @param output_dir (character) directory to save TIF files if write files. Default = NULL
#' @param overwrite (logical) whether to overwrite files in the output directory. Default = FALSE.
#' @param progress_bar (logical) whether to show progress bar. Default = TRUE.
#'
#' @return
#' A list of SpatRasters containing the interpolated bioclimatic variables for
#' each period.
#' @importFrom progress progress_bar
#' @importFrom terra writeRaster
#' @export
#' @details
#' The interpolated bioclimatic variables are constructed using anomalies between
#' bioclimatic layers of the present and the Last Glacial Maximum (LGM), scaled
#' relative to the Plio-Pleistocene global mean temperature curve derived from benthic
#' stable oxygen isotope ratios (see `data("scalings")`, `?scalings` and references
#' for details).The function allows interpolation across 538 scenarios, spanning
#' from the LGM (~20,000 years before present) to 5.4 million years ago.
#'
#' @references
#' Gamisch, A. (2019). Oscillayers: A dataset for the study of climatic
#' oscillations over Plio‐Pleistocene time‐scales at high spatial‐tempora
#' l resolution. Global Ecology and Biogeography, 28(11), 1552-1560.
#'
#'
oscillayer <- function(current_variables, lgm, final_period = 120,
                      precipitation_var = paste0("bio_", 12:19), #To replace negative values by 0
                      write_files = FALSE,
                      output_dir = NULL,
                      overwrite = FALSE,
                      progress_bar = TRUE) {

  #Import scalings
  scalings <- grinnell::scalings$scaling_factor[-1]

  #Get vector of periods (in thousand of years before present)
  periods <- seq(30, final_period, 10)

  #Calculate delta (current - lgm)
  delta <- current_variables - lgm

  if (progress_bar) {
    # Replicates
    total_steps <- length(periods) * terra::nlyr(delta)
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent eta::eta [:elapsed]",
      total = total_steps,
      clear = FALSE
    )
  }
  res_period <- list()
  res_l <- list()
  for(f in 1:length(periods)){
    scaling_i <- scalings[f]
    for(l in names(delta)){
      if (progress_bar) {
        pb$tick()
      }
      res_l[[l]] <- sum(c((delta[[l]] * scaling_i), lgm[[l]]))
    }
    res_period[[f]] <- rast(res_l)

    if(any(terra::minmax(res_period[[f]][[precipitation_var]])[1,] < 0)){
      res_period[[f]][[precipitation_var]] <- (res_period[[f]][[precipitation_var]]
                                               > 0) * res_period[[f]][[precipitation_var]]}
  }
  names(res_period) <- paste0(periods, "kya")

  if(write_files){
    for(i in names(res_period)){
      terra::writeRaster(res_period[[i]],
                         paste0(output_dir, "/", i, ".tif"),
                         overwrite =  overwrite)
    }
  }

  return(res_period)
}

# #Get scaling factors from Scaling_factor_TabS1.csv in https://datadryad.org/stash/dataset/doi:10.5061/dryad.27f8s90
# scalings <- read.csv("c:/Users/wever/Downloads/Input/Scaling_factor_TabS1.csv")
# scalings$SL..m. <- NULL
# colnames(scalings) <- c("Myr_bp",
#                         "Time_period_10ky_bp",
#                         "Surface_temperature",
#                         "delta_LGM",
#                         "scaling_factor")
# usethis::use_data(scalings, overwrite = TRUE)

# #Test function
# current_variables <- rast("../paleoclim_10/cur.tif")
# lgm <- rast("../paleoclim_10/lgm.tif")
# output_dir = "../paleoclim_10/"
# write_files = TRUE
