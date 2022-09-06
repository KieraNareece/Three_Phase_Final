# ======================================================================================================================
# Script Name: Optimisation.py
# File Name: Three_Phase_Final
# Date: 06-09-2022
# Creator: K.N. Lambrecht
# Purpose: This script is used for solution of the system of ODEs and the optimisation of parameters
# ======================================================================================================================

from Main_Script import *
from Function import *
import time

# Initial solution of ODE system with initial parameter set

sol = speed_up(p, x)
solution = vector_to_namespace(sol, p.svar_list)
v = simple_calc(ts, solution, p, time_array)

# This is the optimisation using a least squares method

# A new name space for parameters is generated for the optimisation
gs = p.param_regression
gl = p.param_names_reg
g = (namespace_to_vector(gs, gl))

# This is used as a placeholder for when bootstrapping occurs. It is used to keep the model prediction vector in the
# correct order when bootstrapping does not occur
optimise_index = np.arange(0, len(exp_vector))

# Regression using scipy.optimize.least_squares
start = time.time()
# a = least_squares(objective, g, method='trf', jac='3-point', bounds=bounds, x_scale='jac',
#                   tr_solver='lsmr', args=(p, x, ts, time_array, gl, determinate, exp_vector, stdev_vector,
#                                           optimise_index), verbose=2)
stop = time.time()
print(stop-start)

# Once the optimisation has finished, the new parameters are converted into a namespace and replace the initial guesses
# for the model

new_p = vector_to_namespace(g, p.param_names_reg)
p.param_regression = new_p
sol2 = speed_up(p, x)
solution2 = vector_to_namespace(sol, p.svar_list)
v2 = simple_calc(ts, solution, p, time_array)


subplots(ts, time_array, solution, v, solution2, v2, average, sd)
