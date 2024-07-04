#' Get and predict models to current and past conditions
#'
#' @description get_models uses ellipsoid_suitability to produce a model of
#' environmental suitability based on the Mahalanobis distance to the centroid
#' of environmental values characterized by species occurrences. It then
#' predicts the models for past conditions.
#'
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
#'
#' @return
#' A list containing:
#' - models: suitability models (SpatRasters) for current and past scenarios.
#' - initial_m: a SpatVector representing the initial M.
#' - spatial_points: a SpatVector containing the spatialized records.
#'
#' @details
#' If `initial_m_buffer` is not NULL, the current variables are cropped using an
#' initial M. This initial M is based on a minimum convex polygon around the
#' records, plus a buffer of x kilometers defined by `initial_m_buffer`. If
#' `initial_m_buffer` is NULL, the current variables are not cropped.
#'
#' A principal component analysis (PCA) is performed with `current_variables`.
#' The first three principal components are then used to calculate the
#' suitability layer for dispersal simulations. Suitability values are derived
#' from an ellipsoid envelope model created based on occurrence records and
#' principal components. The ellipsoid model is used because it is a simple yet
#' reliable representation of a species' ecological niche that does not require
#' background or pseudo-absences.
#'
#' Lastly, the suitability models are projected to past conditions.
#'
#' @importFrom terra extract convHull buffer prcomp predict disagg as.polygons
#' trim rast crop
#' @export
#' @usage get_models(data, long, lat, current_variables, variables = NULL,
#'                   scale = TRUE, center = TRUE, projection_dir, periods,
#'                   pattern = ".tif", initial_m_buffer = NULL,
#'                   suitability_threshold = 5)
get_models <- function(data, long, lat,
                       current_variables,
                       variables = NULL,
                       scale = TRUE,
                       center = TRUE,
                       projection_dir,
                       periods,
                       pattern = ".tif",
                       initial_m_buffer = NULL,
                       suitability_threshold = 5) {

  #Subset current variables if necessary
  if(!is.null(variables)){
    current_variables <- current_variables[[variables]]
  }

  #Get occurrences
  xy <- data[,c(long, lat)]
  #Convert points to spatvector
  xys <- vect(xy, geom = c(x = long, y = lat), crs = "+init=epsg:4326")
  #mapview(xys)

  #Remove NAs from PCAs
  remove_na <- which(is.na(terra::extract(current_variables[[1]],
                                          xy, ID = FALSE, raw = TRUE)))
  if(length(remove_na) > 0){
    data <- data[-remova_na,]
    xy <- xy[-remove_na,]}

  #Set initial m
  if(!is.null(initial_m_buffer)){
    #Get minimum convex polygon and add distance
    big_m <- terra::convHull(xys)
    #Add buffer
    big_m <- terra::buffer(big_m, width = initial_m_buffer*1000)
    #Cut to continent
    continent <- as.polygons(!is.na(current_variables[[1]]) * 1)
    continent <- continent[continent$bio_1 == 1]
    big_m <- crop(big_m, continent)
    #Disaggregate big_m
    big_m <- disagg(big_m)

    #Cut variables of the present
    current_variables <- crop(current_variables, big_m, mask = TRUE) }


  #Do PCA
  pca <- terra::prcomp(current_variables, scale = scale, center = center)
  #Predict and get 3 first axes
  p_pca <- terra::predict(current_variables, pca, index = 1:3)

  #Get ellipse model
  m <- ellipsoid_suitability(data = xy,
                             suitability_threshold = suitability_threshold,
                             variables = p_pca)

  #Get model of present
  m_p <- m$suitability_layer

  #Cut big m using model of present
  new_big_m <- terra::disagg(terra::as.polygons(terra::trim(!is.na(m_p) * 1,
                                                            0)))

  past <- list()
  for(i in 1:length(periods)){
    #Get period i
    period_x <- periods[i]
    p_i <- terra::rast(paste0(projection_dir, "/", period_x, pattern))
    #Cut variables
    p_i <- terra::crop(p_i, m_p, mask = TRUE)
    #Project pca
    p_i <- terra::predict(p_i, pca, index = 1:3)
    #Predict model
    m_i <- predict_esuitability(ellipsoid_model = m,
                                variables = p_i,
                                suitability_threshold = suitability_threshold)
    past[[i]] <- m_i$suitability_layer}
    #mbin_i <- (m_i$suitability_layer > thr) * m_i$suitability_layer
    #NAflag(mbin_i) <- 0

  names(past) <- periods

  #Save all results in a list
  return(list(models = rast(c(list(cur = m_p), past)),
              initial_m = new_big_m,
              data = data))
}

