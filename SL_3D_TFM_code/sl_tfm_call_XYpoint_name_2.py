#~~~~~~~~~~~~~ Single-layer TPT-based Traction Force Microscopy ~~~~~~~~~~~~~~~~
#
#Caller script to run the FEniCS solver to compute displacement and stress
#(inc. traction).  Prior to running this script, the 3D single-layer displacements
#should be computed (by TPT) and regularized. The output is saved in .mat files
#prefixed by the data_name_in variable. Inputs (data name, mat props) are set by
# the Matlab script "updatePyRun.m", and should NOT be modified in this script.
#WARNIG: line numbers are hardcoded in the Matlab script.
#
#--- INPUTS --- (set in the file, written via Matlab set up script)
#  data_name_in : name of the input data file containing 3D single-layer displacements
#  data_name_out: name prefix of the output data files
#  load_steps   : number of load steps used to compute the solution
#                  (fewer => faster; greater => more robust convergence)
#  E            : linearized elastic modulus of the substrate gel
#  nu           : Poisson's ratio of the substrate gel
#  thickness    : Estimated gel thickness
#  x_center_norm: normalized x-coord of the center of the deformation
#  y_center_norm: normalized x-coord of the center of the deformation
#
# June, 2019; Alex Landauer
# Franck Lab, Brown University and University of Wisc - Madison
#
from sl_tfm_fenics_dev import sl_tfm_solve #import the solver function

for stepNum in range(2):
	#define inputs for the solver function
	data_name_in = "example_data_XYpoint_name_2_fenicsIn_%03d" % (stepNum)
	data_name_out = "example_data_XYpoint_name_2_fenicsOut_%03d" % (stepNum)
	load_steps = 10
	E = 1500.00000
	nu = 0.45000
	thickness = 70.00000
	x_center_norm = 0.50000
	y_center_norm = 0.50000

	#call the solver function
	sl_tfm_solve(data_name_in,data_name_out,load_steps,E,nu,thickness,x_center_norm,y_center_norm)