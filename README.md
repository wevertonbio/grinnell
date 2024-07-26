An alternative version of grinnell: Dispersal Simulations Based on
Ecological Niches
================
Weverton Trindade

- [Project description](#project-description)
- [What is new in this alternative
  version?](#what-is-new-in-this-alternative-version)
- [Installation](#installation)
- [Workflow](#workflow)
  - [1 Download variables from
    rpaleoclim](#1-download-variables-from-rpaleoclim)

<hr>

## Project description

This is an alternative version of the [grinnell R
package](https://github.com/fmachados/grinnell). Please, read and cite
the [original
paper](https://escholarship.org/content/qt8hq04438/qt8hq04438.pdf).

**grinnell** is an R package to simulate dispersal, colonization, and
accessibility based on niche estimations. One of the main algorithms
implemented here is the simulation of species accessible areas (M),
which can be used as calibration areas in Ecological Niche Models (ENM)
and Species Distribution Models (SDM).

To simulate M, grinnell uses the same inputs needed by several ENM and
SDM: clean occurrences of the study species and a set of climatic layers
(rasters).

## What is new in this alternative version?

While the original version constructs glacial-interglacial climate
conditions based on interpolations between current and LGM climate
conditions, this alternative version allows for the use of climatic
variables representing glacial-interglacial climate conditions provided
externally by other sources, such as
[PaleoClim](https://cran.r-project.org/web/packages/rpaleoclim/vignettes/rpaleoclim.html)
and
[Oscillayers](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12979).

## Installation

To install and call this alternative version of `grinnell` use:

``` r
if (!require("remotes")) {
  install.packages("remotes")
}

remotes::install_github("wevertonbio/grinnell")
library(grinnell)
```

<hr>

<br>

## Workflow

### 1 Download variables from rpaleoclim

First, we need to download the variables representing climatic
conditions between the current period and LGM period from
[PaleoClim](https://cran.r-project.org/web/packages/rpaleoclim/vignettes/rpaleoclim.html).
To do this, install the *rpaleoclim* R package and download the
variables:

``` r
#### Download paleoclim data with rpaleoclim ####
if (!require("rpaleoclim")) {
  install.packages("rpaleoclim")
}
library(rpaleoclim)
library(terra)
library(pbapply)

#Create directory to save layers
dir.create("paleoclim_10")

#Import neotropic vector to cut variables
#You can use any other vector
neot <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/Neotropic.gpkg")
#Add buffer around the polygon
neot <- buffer(neot, width = 100*1000)

#Get times to download
times <- c("cur", #Current (1979 – 2013)
           "lh", #Late Holocene: Meghalayan 4.2-0.3 ka
           "mh", #Mid Holocene: Northgrippian   8.326-4.2 ka
           "eh", #Early Holocene: Greenlandian  11.7-8.326 ka
           "yds", #Pleistocene: Younger Dryas Stadial   12.9-11.7 ka
           "ba",    #Pleistocene: Bølling-Allerød   14.7-12.9 ka
           "hs1", #Pleistocene: Heinrich Stadial 1  17.0-14.7 ka
           "lgm", #Pleistocene: Last Glacial Maximum    ca. 21 ka
           "lig") #Pleistocene: Last Interglacial   ca. 130 ka

#### Download variables ar 10arcmmin of resolution ####
pblapply(times, function(i){
  try({
    p <- paleoclim(period = i,
                   resolution = "10m", #10arc-min of resolution
                   as = "terra", 
                   region = neot)
    #Save variables
    writeRaster(p, paste0("paleoclim_10", i, ".tif"))
  }) #End of try
})
```

Note that there is a significant interval between the LGM (ca. 21ka) and
the LIG (ca. 130 ka). We can fill this gap using the same method
employed by
[Oscillayers](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12979).
Oscillayers constructs environmental variables representing climatic
conditions over the last 5.4 million years at intervals of 10ka. It is
based on interpolated anomalies between bioclimatic layers of the
present and the Last Glacial Maximum (LGM), scaled relative to the
Plio-Pleistocene global mean temperature curve derived from benthic
stable oxygen isotope ratios. This alternative version of grinnell
includes a function to generate these layers. Let’s use the current and
LGM variables downloaded from PaleoClim as inputs for this function:

``` r
#Run this only once to generate the variables between LGM and the desired time.
current_variables <- rast("paleoclim_10/cur.tif")
lgm <- rast("paleoclim_10/lgm.tif")
oscillayer(current_variables = current_variables,
           lgm = lgm,
           final_period = 120, #Final period, here 120 ka
           precipitation_var = c("bio_12", "bio_13", "bio_14", #Identify precipitation variables to make corrections
                                   "bio_15", "bio_16", "bio_17"),
           write_files = TRUE,
           output_dir = "paleoclim_10", #Write the layers in the same folder where we save rpaleoclim variables
           overwrite = TRUE,
           progress_bar = TRUE)
```

Check the folder set in *output_dir* and see the new variables.

Now, we can use these variables to run the simulations:

``` r
#Load current variables
current_variables <- rast("paleoclim_10/cur.tif")
#Set the sequence of periods, from the most recent to the oldest
periods <- c("lh", "mh", "eh", "yds", "ba", "hs1", "lgm", "30kya", "40kya",
             "50kya", "60kya", "70kya", "80kya", "90kya",
             "100kya", "110kya", "120kya","lig")
#Set the bioclimatic variables to be used
variables <- c("bio_1", "bio_10", "bio_11", "bio_12", "bio_13", "bio_14",
               "bio_15", "bio_16", "bio_17","bio_2", "bio_3",
               "bio_4", "bio_5", "bio_6", "bio_7")
#Set the folder with the bioclimatic variables
projection_dir <- "paleoclim_10"

#Grid of combinations to test M: dispersal events and kernel_spread
#See ?m_simulations for more details
comb_grid <- data.frame(dispersal_events = c(5, 10, 15, 20, 25, 30, 35),
                        kernel_spread = c(2, 2, 2, 2, 3, 4, 6))

#Run M simulation
m <- m_simulations(data = occ, #Dataframe with longitude and latitude of the records
                   long = "x", lat = "y", #Column names with longitude and latitude
                   current_variables = current_variables,
                     variables = variables,
                     scale = TRUE,
                     center = TRUE,
                     projection_dir = projection_dir,
                     periods = periods,
                     pattern = ".tif",
                     initial_m_buffer = 500,
                     suitability_threshold = 5,
                     starting_proportion = 0.75,
                     proportion_to_disperse = 1,
                     sampling_rule = "random",
                     dispersal_kernel = "normal",
                     kernel_spread = 2,
                     max_dispersers = 4,
                     dispersal_events = 25,
                     comb_grid = comb_grid,
                     replicates = 3,
                     threshold = 5,
                     set_seed = 1,
                     skip_extinction = TRUE,
                     results_by_event = FALSE,
                     results_by_scenario = FALSE,
                     remove_m_without_records = TRUE,
                     extra_buffer = 50,
                     progress_bar = TRUE,
                     verbose = TRUE)
```

You can use the information in **m\$summary** to select the optimal
combination of M, based on the number of disjoint polygons and the
number of records outside the final M.

If you set **results_by_event** or **results_by_scenario** to TRUE, you
can use these results to make a GIF showing the dispersal simulation by
event or by scenario:

``` r
gif_dispersion(m$m_by_scen, gif_file = "Dispersion_by_scenario.gif",
               width = 2200,
               height = 1500,
               res = 300,
               delay = 0.1,
               loop = FALSE,
               progress_bar = TRUE,
               verbose = TRUE)
gif_dispersion(m$m_by_event, gif_file = "Dispersion_by_event.gif",
               width = 2200,
               height = 1500,
               res = 300,
               delay = 0.1,
               loop = FALSE,
               progress_bar = TRUE,
               verbose = TRUE
               
```
