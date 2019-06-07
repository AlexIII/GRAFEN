Specify file for the output field:
-dat2D dat2DfiledExample.dat		- XYV file, where V is the gravity filed value. XY - arbitrary points in Gauss-Kruger projection (km). V is ignored as input. See example dat2DfiledExample.dat
or
-dat3D dat3DfiledExample.dat		- XYZV file, where V is the gravity filed value. XY - arbitrary points in Gauss-Kruger projection (km). Z and V are ignored as input. Z will be set to H (see bellow). See example dat2DfiledExample.dat
or
-grd7 grdField.grd			- Grid dementions are used as file input. All coorditates in Gauss-Kruger projection (km). Old grid valuse are ignored and will be rewritten.

Specify input parameters:
-Hf 0.00001				- Height over the Ellipsoid on which the field is being calculated (km). (It is not recommended to pass exactly 0)
-dens dens_model			- directory with layerd density model (set of grd7 files)
-Hfrom -81 				- Depth of the lower grid layer (density model)
-Hto 0 					- Depth of the upper grid layer (density model)
-Hn 81 					- number of layers of the density model (must be same as amount of files in dens_model)
-l0 57 					- central meridean for Gauss-Kruger projection
-DPR 180 				- radius of point-field replacement (in km)
-toRel					- convert input density model to relative values
-noInvFileOrder				- don't invert the file order of density model
-transposeSolver			- solve gravity problem with transposed forward gravity field operator. Files in dens_model will be rewritten, "output field" is now treeted as input.
