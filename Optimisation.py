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

sol = speed_up(p, x1)
solution = vector_to_namespace(sol, p.svar_list)
v = simple_calc(ts, solution, p, time_array)

# This is the optimisation using a least squares method

# This is used as a placeholder for when bootstrapping occurs. It is used to keep the model prediction vector in the
# correct order when bootstrapping does not occur
optimise_index_1 = np.arange(0, len(exp_vector1))

# Regression using scipy.optimize.least_squares
start = time.time()
a1 = least_squares(objective, g, bounds=bounds, jac='3-point', x_scale='jac', loss='soft_l1', tr_solver='lsmr',
                   args=(p, x1, ts, time_array, gl, determinate1, exp_vector1, stdev_vector1, optimise_index_1),
                   verbose=2)

stop = time.time()
print(a1, stop-start)

optimise_index_2 = np.arange(0, len(exp_vector2))

# Regression using scipy.optimize.least_squares
start = time.time()
a2 = least_squares(objective, g, bounds=bounds, jac='3-point', x_scale='jac', loss='soft_l1', tr_solver='lsmr',
                   args=(p, x2, ts, time_array, gl, determinate2, exp_vector2, stdev_vector2, optimise_index_2),
                   verbose=2)

stop = time.time()
print(a2, stop-start)

optimise_index_3 = np.arange(0, len(exp_vector3))

# Regression using scipy.optimize.least_squares
start = time.time()
a3 = least_squares(objective, g, bounds=bounds, jac='3-point', x_scale='jac', loss='soft_l1', tr_solver='lsmr',
                   args=(p, x3, ts, time_array, gl, determinate3, exp_vector3, stdev_vector3, optimise_index_3),
                   verbose=2)

stop = time.time()
print(a3, stop-start)

optimise_index_4 = np.arange(0, len(exp_vector4))

# Regression using scipy.optimize.least_squares
start = time.time()
a4 = least_squares(objective, g, bounds=bounds, jac='3-point', x_scale='jac', loss='soft_l1', tr_solver='lsmr',
                   args=(p, x4, ts, time_array, gl, determinate4, exp_vector4, stdev_vector4, optimise_index_4),
                   verbose=2)

stop = time.time()
print(a4, stop-start)

optimise_index_5 = np.arange(0, len(exp_vector5))

# Regression using scipy.optimize.least_squares
start = time.time()
a5 = least_squares(objective, g, bounds=bounds, jac='3-point', x_scale='jac', loss='soft_l1', tr_solver='lsmr',
                   args=(p, x5, ts, time_array, gl, determinate5, exp_vector5, stdev_vector5, optimise_index_5),
                   verbose=2)

stop = time.time()
print(a5, stop-start)

optimise_index_6 = np.arange(0, len(exp_vector6))

# Regression using scipy.optimize.least_squares
start = time.time()
a6 = least_squares(objective, g, bounds=bounds, jac='3-point', x_scale='jac', loss='soft_l1', tr_solver='lsmr',
                   args=(p, x6, ts, time_array, gl, determinate6, exp_vector6, stdev_vector6, optimise_index_6),
                   verbose=2)

stop = time.time()
print(a6, stop-start)

optimise_index_7 = np.arange(0, len(exp_vector7))

# Regression using scipy.optimize.least_squares
start = time.time()
a7 = least_squares(objective, g, bounds=bounds, jac='3-point', x_scale='jac', loss='soft_l1', tr_solver='lsmr',
                   args=(p, x7, ts, time_array, gl, determinate7, exp_vector7, stdev_vector7, optimise_index_7),
                   verbose=2)

stop = time.time()
print(a7, stop-start)

optimise_index_8 = np.arange(0, len(exp_vector8))

# Regression using scipy.optimize.least_squares
start = time.time()
a8 = least_squares(objective, g, bounds=bounds, jac='3-point', x_scale='jac', loss='soft_l1', tr_solver='lsmr',
                   args=(p, x8, ts, time_array, gl, determinate8, exp_vector8, stdev_vector8, optimise_index_8),
                   verbose=2)

stop = time.time()
print(a8, stop-start)

optimise_index_9 = np.arange(0, len(exp_vector9))

# Regression using scipy.optimize.least_squares
start = time.time()
a9 = least_squares(objective, g, bounds=bounds, jac='3-point', x_scale='jac', loss='soft_l1', tr_solver='lsmr',
                   args=(p, x9, ts, time_array, gl, determinate9, exp_vector9, stdev_vector9, optimise_index_9),
                   verbose=2)

stop = time.time()
print(a9, stop-start)

# Once the optimisation has finished, the new parameters are converted into a namespace and replace the initial guesses
# for the model

new_p = vector_to_namespace(a1.x, p.param_names_reg)
p.param_regression = new_p
sol2 = speed_up(p, x1)
solution2 = vector_to_namespace(sol, p.svar_list)
v2 = simple_calc(ts, solution, p, time_array)


subplots(ts, time_array, solution, v, solution2, v2, average, sd)
