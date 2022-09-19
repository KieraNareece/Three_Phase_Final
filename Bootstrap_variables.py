# ======================================================================================================================
# Script Name: Bootstrap.py
# File Name: Three_Phase_Final
# Date: 13-09-2022
# Creator: K.N. Lambrecht
# Purpose: This script is used to perform the likelihood for the parameter sets
# ======================================================================================================================

from Main_Script import *
from Function import *
import numpy as np

num_iterations = 2
bootstrap_indices = np.zeros((num_iterations, len(exp_vector1)), dtype=int)

for i in range(0, num_iterations):
    new_exp, new_sd = [], []
    bstp_index = bill(determinate1, time_array)
    for j in range(0, len(bstp_index)):
        new_exp.append(exp_vector1[bstp_index[j]]), new_sd.append(stdev_vector1[bstp_index[j]])

    bootstrap_indices[i, :] = bstp_index













