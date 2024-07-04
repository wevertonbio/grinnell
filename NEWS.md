# grinnell 0.0.23 - Weverton version

* M_simulationR and dispersal_simulationR: Allowed parallelization of replicates and fixed a bug related to overwriting.
* dispersal_helpers:
	- In nd_sval, replaced apply with rowsums to improve performance.
	- In set_pop, initial_colonized, and suitable_cells, added arguments to identify longitude and latitude.
	- Also fixed a bug when rule = "suitability".
* pca_raster: Fixed bug when var_points has NAs.
* short_helpers: 
	- created new raster using terra::init to improve performance.
	- Fixed a bug with overwriting.
	- Added function make_plot() to help make gifs
* oscillayer and data/scallings: Data and code to interpolate variables from LGM to earlier periods.
* get_models: new function to get and predict models to current and past conditions.
* scenario_wise_simulation_b: alternative version of scenario_wise_simulation to simulate species dispersal processes.
* m_simulations: alternative version of M_simulationR to simulate species dispersal processes.
* gif_dispersion: new function to render animated gifs showing dispersion by event or scenario

# grinnell 0.0.22

* Initial package on github.
