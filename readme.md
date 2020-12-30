# GRAFEN - TOPO

## Description

This is a special version of [GRAFEN](https://github.com/AlexIII/GRAFEN) that allows for gravity field calculation for a single-layer ellipsoidal (or "falt") density model bounded by the topography from the above (height map). For now, only constant density distribution is supported. 
This version is a work in progress.

For the build instructions please refer to the `master` branch `readme.md`.


## Program arguments
Specify file for the output field:

`-gravGrd7 grdField.grd`			- Surfer grd7 file. Grid dimensions are used as file input. All coordinates in Gauss-Kruger projection (km). Old grid values are ignored and will be rewritten.

Specify input parameters:
`-topoGrd7 grdTopo.grd`			- Surfer grd7 file. Height map (upper boundry of the model). Grid dimensions are in Gauss-Kruger projection (km).

`-dens VAL`			- Model constant density in g/cm^3. A number.

`-Req VAL`			- (optional) Equatorial radius, km. (6378.245km if not specified)
`-Rpol VAL`			- (optional) Polar radius, km. (6356.863km if not specified)

`-DPR 180` 			-(optional) Radius of point-field replacement (in km). If you don't specify this option, replacement radius will be automatically deduced, based on condition that the output field accuracy won't be reduced more than by 0.1%. Set to '-1' to disable (default).

`-nx`			      - (optional) Calculate gravity field along x normal.
`-ny`			      - (optional) Calculate gravity field along y normal.
`-nz`			      - (optional) Calculate gravity field along z normal.

`-flat` 				- Calculate gravity field for 'flat' model.


## License

This software is distributed under MIT License. © Alexander Chernoskutov, Denis Byzov


## Acknowledgements

Development of this software was supported by the Russian Foundation for Basic Research (project no. 20-05-00230_a). 

