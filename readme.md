# GRAFEN - TOPO

## Description

This is a special version of [GRAFEN](https://github.com/AlexIII/GRAFEN) that allows for gravity field calculation for a single-layer ellipsoidal (or "falt") density model bounded by the topography from the above (height map). For now, only constant density distribution is supported. 
This version is a work in progress.

For the build instructions please refer to the `master` branch `readme.md`.


## Program arguments

Specify data files

`-grd7 grdField.grd`		- Gravity field field input-output. Surfer grd7 file. Grid dimensions are used as file input. All coordinates in Gauss-Kruger projection (km). Old grid values are ignored and will be rewritten.

`-dens dens_model`			- Directory with layered density model (set of grd7 files). The top layer treated as density distribution in the topography layer.

`-topoHeightGrd7 topo.grd`  - Surfer grd7 file. Height map (upper boundry of the model) in km. Grid dimensions are in Gauss-Kruger projection (km).

Specify input parameters:

`-fieldOnTopo`				- (optional) Calculate field on the topography heights.

`-Hf 0.00001`				- Height over the Ellipsoid on which the field is being calculated (km). This parameter is ignored if `-dat3D` was specified. It is not recommended to pass exactly 0.

`-Hfrom -81` 				- Depth of the lower grid layer (density model)

`-Hto 0` 					- Depth of the upper grid layer (density model)

`-Hn 81` 					- Number of layers of the density model (must be same as amount of files in dens_model)

`-l0 57` 					- Central meridian for Gauss-Kruger projection

`-DPR 180` 				    - (optional) Radius of point-field replacement (in km). If you don't specify this option, replacement radius will be automatically deduced, based on condition that the output field accuracy won't be reduced more than by 0.1%.

`-toRel`					- (optional) Convert input density model to relative values

`-noInvFileOrder`			- (optional) Don't invert the file order of density model

`-transposeSolver`		    - (optional) Solve gravity problem with transposed forward gravity field operator. Files in dens_model will be rewritten, "output field" is now treated as input.

`-Req VAL`			        - (optional) Equatorial radius, km. (6378.245km if not specified)

`-Rpol VAL`			        - (optional) Polar radius, km. (6356.863km if not specified)


## License

This software is distributed under MIT License. © Alexander Chernoskutov, Denis Byzov


## Acknowledgements

Development of this software was supported by the Russian Foundation for Basic Research (project no. 20-05-00230_a). 

