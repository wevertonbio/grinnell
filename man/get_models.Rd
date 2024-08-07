% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_models.R
\name{get_models}
\alias{get_models}
\title{Get and predict models to current and past conditions}
\usage{
get_models(data, long, lat, current_variables, variables = NULL,
                  scale = TRUE, center = TRUE, projection_dir, periods,
                  pattern = ".tif", initial_m_buffer = NULL,
                  suitability_threshold = 5)
}
\arguments{
\item{data}{(data.frame) occurrence records for the species of interest to be
used to run the simulation. Columns must include species name, longitude,
and latitude}

\item{long}{(character) column name in data with longitude}

\item{lat}{(character) column name in data with latitude}

\item{current_variables}{(SpatRaster) environmental variables representing
"current" conditions.}

\item{variables}{(character) variable names to subset current variables. If
NULL (default), use all variables.}

\item{scale}{(logical) whether or not to scale variables while performing
principal component analyses.
Default = TRUE.}

\item{center}{(logical) whether or not to center variables while performing
principal component analyses.
Default = TRUE.}

\item{projection_dir}{(character) directory with raster files representing
past conditions.}

\item{periods}{(character) names of the periods representing past conditions,
order from the most recent to oldest.}

\item{pattern}{(character) format of the raster files representing past
conditions. Default = ".tif".}

\item{initial_m_buffer}{(numeric) width (in kilometers) of the buffer to be
added to the initial M (See details). If NULL (default), the initial M will
have exactly the same extent as the current variables.}

\item{suitability_threshold}{(numeric) value (in percentage) to be used as
threshold for suitability. Below this value environments are considered
unsuitable. Default = 5.}
}
\value{
A list containing:
\itemize{
\item models: suitability models (SpatRasters) for current and past scenarios.
\item initial_m: a SpatVector representing the initial M.
\item spatial_points: a SpatVector containing the spatialized records.
}
}
\description{
get_models uses ellipsoid_suitability to produce a model of
environmental suitability based on the Mahalanobis distance to the centroid
of environmental values characterized by species occurrences. It then
predicts the models for past conditions.
}
\details{
If \code{initial_m_buffer} is not NULL, the current variables are cropped using an
initial M. This initial M is based on a minimum convex polygon around the
records, plus a buffer of x kilometers defined by \code{initial_m_buffer}. If
\code{initial_m_buffer} is NULL, the current variables are not cropped.

A principal component analysis (PCA) is performed with \code{current_variables}.
The first three principal components are then used to calculate the
suitability layer for dispersal simulations. Suitability values are derived
from an ellipsoid envelope model created based on occurrence records and
principal components. The ellipsoid model is used because it is a simple yet
reliable representation of a species' ecological niche that does not require
background or pseudo-absences.

Lastly, the suitability models are projected to past conditions.
}
