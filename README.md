# consopt
Comm Satellite Constellation design tool developed for my Master's dissertation.

This tool couples the [MGEO](https://github.com/ronisbr/MGEO.jl/) algorithm with [SatelliteToolbox](https://github.com/JuliaSpace/SatelliteToolbox.jl) to obtain a Pareto Frontier for the problem of maximum time gap
for a set of data collection stations, receiving stations and constellation.

Main file is run_opt.jl, used to call the optimization loop.
functions_opt.jl includes data initialization and functions used overall.
standalonerun file is used to check overall calculations apart from the optimization loop.
