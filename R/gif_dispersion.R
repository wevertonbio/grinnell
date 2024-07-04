#' Render animated gifs showing dispersion by event or scenario
#'
#' @param rasters (SpatRasters) Combined layers representing the accessed, colonized,
#' and suitable-but-inaccessible cells by scenario or event. This is the output of
#' `m_simulations` or `scenario_wise_simulation_b` when `results_by_event` or
#' `results_by_scenario` are set to TRUE.
#' @param gif_file (character) Path and name of the GIF file to be created.
#' @param width (numeric) gif width in pixels. Default = 2200.
#' @param height (numeric) gif height in pixels. Default = 1500.
#' @param res (numeric) gif resolution in pixels. Default = 300.
#' @param delay (numeric) time to show each image in seconds. Default = 0.1.
#' @param loop (logical/numeric) if the gif should be repeated. Set to FALSE
#' (default) to only play once, or a number to indicate how many times to repeat
#' after the first.
#' @param progress_bar (logical) whether to show progress bar. Default = TRUE.
#' @param verbose (logical) whether to display messages.
#'
#' @export
#' @importFrom gifski save_gif
#' @importFrom terra plot nlyr
#' @importFrom progress progress_bar
gif_dispersion <- function(rasters, gif_file, width = 2200,
                           height = 1500,
                           res = 300,
                           delay = 0.1,
                           loop = FALSE,
                           progress_bar = TRUE,
                           verbose = TRUE){
  gifski::save_gif(make_plot(rasters, progress = progress_bar, verbose),
                   gif_file = gif_file, width = width,
                   height = height ,
                   res = res,
                   delay = delay, loop = loop)
}
