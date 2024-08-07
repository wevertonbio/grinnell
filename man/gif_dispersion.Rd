% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gif_dispersion.R
\name{gif_dispersion}
\alias{gif_dispersion}
\title{Render animated gifs showing dispersion by event or scenario}
\usage{
gif_dispersion(
  rasters,
  gif_file,
  width = 2200,
  height = 1500,
  res = 300,
  delay = 0.1,
  loop = FALSE,
  progress_bar = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{rasters}{(SpatRasters) Combined layers representing the accessed, colonized,
and suitable-but-inaccessible cells by scenario or event. This is the output of
\code{m_simulations} or \code{scenario_wise_simulation_b} when \code{results_by_event} or
\code{results_by_scenario} are set to TRUE.}

\item{gif_file}{(character) Path and name of the GIF file to be created.}

\item{width}{(numeric) gif width in pixels. Default = 2200.}

\item{height}{(numeric) gif height in pixels. Default = 1500.}

\item{res}{(numeric) gif resolution in pixels. Default = 300.}

\item{delay}{(numeric) time to show each image in seconds. Default = 0.1.}

\item{loop}{(logical/numeric) if the gif should be repeated. Set to FALSE
(default) to only play once, or a number to indicate how many times to repeat
after the first.}

\item{progress_bar}{(logical) whether to show progress bar. Default = TRUE.}

\item{verbose}{(logical) whether to display messages.}
}
\description{
Render animated gifs showing dispersion by event or scenario
}
