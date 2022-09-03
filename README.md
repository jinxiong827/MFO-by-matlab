# MFO-by-matlab

Multi-fan optimization framework by Matlab and ANSYS CFX,
"Rapid, concurrent and optimal method for the design of multi-fan using CFD and surrogate-assisted optimization"

10 fans optimization is taken as example.

/MFO/initial_generate.m is the initialization part 

/MFO/opt_start.m is the main loop

Other functions:
/MFO/blade_profile.m, arcpoint.m, included_angle.m and bezier.m: calculate the coordinates of blade points

/MFO/cfxpre.m: write the script for CFD boundary condition setting

/MFO/turbogrid.m: write  the script for generating grid

/MFO/cfxpost.m: write the script for CFD post-processing

/MFO/load_data.m: load CFD data

/MFO/median_value.m: Eq. 4 in paper

/RVEA is refered to PlatEMO(https://github.com/BIMK/PlatEMO)

Kriging toolbox is DACE(http://www2.imm.dtu.dk/projects/dace/).
